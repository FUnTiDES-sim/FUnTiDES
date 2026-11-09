"""Regular rectangular Gmsh blocks with conforming local quad transitions."""
import argparse
import csv
import json
import math
from pathlib import Path

import gmsh
import numpy as np

# Results and reports of the validations
RESULTS = Path(__file__).resolve().parents[3] / "documentation" / "acoustic" / "local-time-step"


def generate(config, output):
    domain, inc, mesh = config["domain"], config["inclusion"], config["mesh"]
    lx, lz, interface = (float(domain[k]) for k in ("lx_m", "lz_m", "interface_z_m"))
    patch = mesh["inclusion_patch"]
    x0, x1, z0 = (float(patch[k]) for k in ("x_min", "x_max", "z_min"))
    z1 = float(mesh["interface_transition_bottom_z_m"])
    if not (0 < x0 < inc["x_min"] < inc["x_max"] < x1 < lx and
            0 < z0 < inc["z_min"] < inc["z_max"] < z1 < interface < lz):
        raise ValueError("The inclusion must lie inside its patch, below the interface transition")
    if not (mesh["design_fmax_hz"] > 0 and mesh["elements_per_wavelength"] > 0 and mesh["transition_m"] > 0):
        raise ValueError("Frequency, resolution and transition length must be positive")
    materials = [(name, config[name]) for name in ("upper_layer", "lower_layer", "inclusion")]
    if any(mat["vp_m_s"] <= 0 or mat["rho"] <= 0 for _, mat in materials):
        raise ValueError("Velocity and density must be positive")
    if len({mat["physical_tag"] for _, mat in materials}) != 3:
        raise ValueError("Material physical tags must be distinct")
    targets = {name: mat["vp_m_s"] / (mesh["design_fmax_hz"] * mesh["elements_per_wavelength"])
               for name, mat in materials}
    # All-quad subdivision halves each edge, also on transfinite quads.
    coarse = 2 * targets["lower_layer"]
    if not all(np.isclose(v / coarse, round(v / coarse)) for v in (x0, x1, lx, z0, z1)):
        raise ValueError("Coarse block coordinates must be multiples of twice the lower-layer target size")
    output.mkdir(parents=True, exist_ok=True)
    gmsh.initialize(["gmsh", "-nopopup"])
    try:
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.model.add("structured_velocity_adapted_bilayer")
        geo = gmsh.model.geo
        points, lines, curves, blocks = {}, {}, {}, []

        def point(xz):
            key = tuple(xz)
            if key not in points:
                points[key] = geo.addPoint(key[0], 0, key[1])
            return points[key]

        def line(a, b):
            a, b = tuple(a), tuple(b)
            if (a, b) in lines:
                return lines[a, b]
            if (b, a) in lines:
                return -lines[b, a]
            tag = geo.addLine(point(a), point(b))
            lines[a, b] = tag
            curves[tag] = (a, b)
            return tag

        def loop(coords):
            return [line(a, b) for a, b in zip(coords, coords[1:] + coords[:1])]

        def surface(name, material, coords, *, holes=(), corners=None, target=None):
            edges = loop(coords)
            loops = [geo.addCurveLoop(edges)]
            for hole in holes:
                loops.append(geo.addCurveLoop(loop(list(reversed(hole)))))
            tag = geo.addPlaneSurface(loops)
            blocks.append(dict(name=name, material=material, surface=tag, edges=edges,
                               corners=[point(v) for v in (corners or coords)], target=target))

        def rect(xa, xb, za, zb):
            return [(xa, za), (xb, za), (xb, zb), (xa, zb)]

        surface("upper_regular", "upper_layer", rect(0, lx, interface, lz), target=targets["upper_layer"])
        surface("lower_left_regular", "lower_layer",
                [(0, 0), (x0, 0), (x0, z0), (x0, z1), (0, z1)],
                corners=rect(0, x0, 0, z1), target=targets["lower_layer"])
        surface("lower_right_regular", "lower_layer",
                [(x1, 0), (lx, 0), (lx, z1), (x1, z1), (x1, z0)],
                corners=rect(x1, lx, 0, z1), target=targets["lower_layer"])
        surface("lower_bottom_regular", "lower_layer", rect(x0, x1, 0, z0), target=targets["lower_layer"])
        inclusion = rect(inc["x_min"], inc["x_max"], inc["z_min"], inc["z_max"])
        surface("inclusion_transition", "lower_layer", rect(x0, x1, z0, z1), holes=[inclusion])
        surface("interface_transition", "lower_layer",
                [(0, z1), (x0, z1), (x1, z1), (lx, z1), (lx, interface), (0, interface)])
        refinement = mesh.get("interior_refinement")
        if refinement:
            xa, xb, za, zb, h = (refinement[k] for k in ("x_min", "x_max", "z_min", "z_max", "element_size_m"))
            if not (inc["x_min"] < xa < xb < inc["x_max"] and inc["z_min"] < za < zb < inc["z_max"] and h > 0):
                raise ValueError("Fine patch must be strictly inside the inclusion")
            inner = rect(xa, xb, za, zb)
            surface("inclusion_inner_transition", "inclusion", inclusion, holes=[inner])
            # Keep exactly the original 40 m boundary discretization. Refinement
            # and its grading live entirely inside the slow material.
            surface("inclusion_fine_regular", "inclusion", inner, target=h)
        else:
            surface("inclusion_regular", "inclusion", inclusion, target=targets["inclusion"])
        geo.synchronize()
        divisions = {}
        if refinement:
            for edge in loop(inclusion):
                a, b = curves[abs(edge)]
                count = math.ceil(math.dist(a, b) / (2 * targets["inclusion"]))
                divisions[abs(edge)] = count
                gmsh.model.mesh.setTransfiniteCurve(abs(edge), count + 1)
        for block in blocks:
            if block["target"] is None:
                continue
            for edge in block["edges"]:
                a, b = curves[abs(edge)]
                count = max(1, math.ceil(math.dist(a, b) / (2 * block["target"]) - 1e-10))
                if abs(edge) in divisions and divisions[abs(edge)] != count:
                    raise ValueError("Incompatible subdivisions on a shared structured edge")
                divisions[abs(edge)] = count
                gmsh.model.mesh.setTransfiniteCurve(abs(edge), count + 1)
            gmsh.model.mesh.setTransfiniteSurface(block["surface"], cornerTags=block["corners"])
            gmsh.model.mesh.setRecombine(2, block["surface"])
        for name, mat in materials:
            tags = [b["surface"] for b in blocks if b["material"] == name]
            gmsh.model.addPhysicalGroup(2, tags, int(mat["physical_tag"]))
            gmsh.model.setPhysicalName(2, int(mat["physical_tag"]), name)
        free, absorbing = [], []
        for tag, (a, b) in curves.items():
            if a[1] == b[1] == lz:
                free.append(tag)
            elif a[1] == b[1] == 0 or a[0] == b[0] == 0 or a[0] == b[0] == lx:
                absorbing.append(tag)
        gmsh.model.addPhysicalGroup(1, free, 1)
        gmsh.model.setPhysicalName(1, 1, "free_surface")
        gmsh.model.addPhysicalGroup(1, absorbing, 2)
        gmsh.model.setPhysicalName(1, 2, "absorbing_sides_bottom")
        # Only the collar and interface strip are unstructured before subdivision.
        # Shared geometric curves enforce identical nodes at every block boundary.
        field = gmsh.model.mesh.field.add("Box")
        for key, value in dict(XMin=inc["x_min"], XMax=inc["x_max"], YMin=-1, YMax=1,
                               ZMin=inc["z_min"], ZMax=inc["z_max"], VIn=2 * targets["inclusion"],
                               VOut=coarse, Thickness=mesh["transition_m"]).items():
            gmsh.model.mesh.field.setNumber(field, key, value)
        gmsh.model.mesh.field.setAsBackgroundMesh(field)
        if refinement:
            # Restrict the fine field to inclusion surfaces, preserving all
            # exterior surface meshes and curve nodes from the baseline.
            fine = gmsh.model.mesh.field.add("Box")
            for key, value in dict(XMin=xa, XMax=xb, YMin=-1, YMax=1, ZMin=za, ZMax=zb,
                                   VIn=2*h, VOut=2*targets["inclusion"],
                                   Thickness=refinement["transition_m"]).items():
                gmsh.model.mesh.field.setNumber(fine, key, value)
            restricted = gmsh.model.mesh.field.add("Restrict")
            gmsh.model.mesh.field.setNumber(restricted, "InField", fine)
            gmsh.model.mesh.field.setNumbers(restricted, "SurfacesList",
                [b["surface"] for b in blocks if b["material"] == "inclusion"])
            combined = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(combined, "FieldsList", [field, restricted])
            gmsh.model.mesh.field.setAsBackgroundMesh(combined)
        for key, value in {"Mesh.MeshSizeFromPoints": 0, "Mesh.MeshSizeFromCurvature": 0,
                           "Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.Algorithm": 6,
                           "Mesh.SubdivisionAlgorithm": 1, "Mesh.ElementOrder": 1,
                           "Mesh.MshFileVersion": 2.2, "Mesh.Binary": 0, "Mesh.SaveAll": 0}.items():
            gmsh.option.setNumber(key, value)
        gmsh.model.mesh.generate(2)
        types, tags, _ = gmsh.model.mesh.getElements(2)
        if list(types) != [3]:
            raise RuntimeError(f"Expected only Q4 elements, got {types}")
        quality = np.asarray(gmsh.model.mesh.getElementQualities(tags[0]))
        if not np.isfinite(quality).all() or quality.min() <= 0:
            raise RuntimeError("Invalid quadrilateral quality")
        counts = {name: 0 for name, _ in materials}
        details = []
        for block in blocks:
            _, element_tags, _ = gmsh.model.mesh.getElements(2, block["surface"])
            count = sum(len(t) for t in element_tags)
            counts[block["material"]] += count
            details.append(dict(name=block["name"], material=block["material"], elements=count,
                                structured=block["target"] is not None))
        gmsh.write(str(output / "bilayer.msh"))
        gmsh.write(str(output / "bilayer.geo_unrolled"))
        config["generation"] = dict(gmsh_version=gmsh.__version__, format="MSH 2.2 ASCII, Q4, XZ plane",
                                    target_final_element_size_m=targets, element_counts=counts, blocks=details,
                                    minimum_quality_minSICN=float(quality.min()),
                                    method="transfinite rectangular blocks; local unstructured transitions; all-quad subdivision")
    finally:
        gmsh.finalize()
    (output / "model.json").write_text(json.dumps(config, indent=2) + "\n")
    with (output / "materials.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["physical_tag", "vp", "rho", "name"])
        for name, mat in materials:
            writer.writerow([mat["physical_tag"], mat["vp_m_s"], mat["rho"], name])
    sim = config["simulation"]
    (output / "simulation.txt").write_text(
        f"{sim['duration_s']} {sim['source_x_m']} {lz-sim['source_depth_m']} {sim['source_sigma_m']} "
        f"{sim['receiver_depth_m']} {sim['cfl']}\n")
    print(json.dumps(config["generation"], indent=2))


if __name__ == "__main__":
    here = Path(__file__).parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=here / "gmsh_bilayer.json")
    parser.add_argument("--output", type=Path, default=RESULTS / "gmsh-structured-results")
    args = parser.parse_args()
    generate(json.loads(args.config.read_text()), args.output)
