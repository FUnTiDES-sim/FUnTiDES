"""Unstructured Q4 meshes of a homogeneous medium with a refined square patch (element size
h outside, h/2 inside), one directory per patch size, for the fine-fraction study of the
local time step (plan step 2). Each directory holds mesh.msh (MSH 2.2, XZ plane: surface 10
= medium, line 1 = free surface on top, line 2 = absorbing sides and bottom), materials.csv
and case.json."""
import argparse
import csv
import json
from pathlib import Path

import gmsh
import numpy as np

# Results and reports of the validations
RESULTS = Path(__file__).resolve().parents[3] / "documentation" / "acoustic" / "local-time-step"


def generate(output, side, h, lx=6000.0, lz=3000.0, vp=2000.0, thickness=200.0, center=(3000.0, 1500.0)):
    output.mkdir(parents=True, exist_ok=True)
    gmsh.initialize(["gmsh", "-nopopup"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.model.add("fine_fraction")
        geo = gmsh.model.geo
        points = [geo.addPoint(x, 0, z, h) for x, z in [(0, 0), (lx, 0), (lx, lz), (0, lz)]]
        edges = [geo.addLine(points[i], points[(i + 1) % 4]) for i in range(4)]   # bottom, right, top, left
        surface = geo.addPlaneSurface([geo.addCurveLoop(edges)])
        geo.mesh.setRecombine(2, surface)
        geo.synchronize()
        if side > 0:
            box = gmsh.model.mesh.field.add("Box")
            cx, cz = center
            for key, value in dict(XMin=cx - side / 2, XMax=cx + side / 2, YMin=-1, YMax=1,
                                   ZMin=cz - side / 2, ZMax=cz + side / 2, VIn=h / 2, VOut=h,
                                   Thickness=thickness).items():
                gmsh.model.mesh.field.setNumber(box, key, value)
            gmsh.model.mesh.field.setAsBackgroundMesh(box)
        for key, value in {"Mesh.MeshSizeFromPoints": 0 if side > 0 else 1, "Mesh.MeshSizeFromCurvature": 0,
                           "Mesh.MeshSizeExtendFromBoundary": 0 if side > 0 else 1,
                           "Mesh.Algorithm": 8, "Mesh.RecombinationAlgorithm": 1,
                           "Mesh.MshFileVersion": 2.2, "Mesh.Binary": 0, "Mesh.SaveAll": 0}.items():
            gmsh.option.setNumber(key, value)
        gmsh.model.addPhysicalGroup(2, [surface], 10, name="medium")
        gmsh.model.addPhysicalGroup(1, [edges[2]], 1, name="free_surface")
        gmsh.model.addPhysicalGroup(1, [edges[0], edges[1], edges[3]], 2, name="absorbing_sides_bottom")
        gmsh.model.mesh.generate(2)
        types, tags, _ = gmsh.model.mesh.getElements(2)
        if list(types) != [3]:
            raise RuntimeError(f"Expected only Q4 elements, got types {list(types)}")
        quality = np.asarray(gmsh.model.mesh.getElementQualities(tags[0]))
        if not np.isfinite(quality).all() or quality.min() <= 0:
            raise RuntimeError("Invalid quadrilateral quality")
        gmsh.write(str(output / "mesh.msh"))
        info = dict(patch_side_m=side, element_size_m=h, patch_element_size_m=h / 2, transition_m=thickness,
                    domain_m=[lx, lz], vp_m_s=vp, elements=int(len(tags[0])),
                    minimum_quality_minSICN=float(quality.min()), gmsh_version=gmsh.__version__,
                    method="frontal-Delaunay for quads + blossom recombination, Box size field")
    finally:
        gmsh.finalize()
    with (output / "materials.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["physical_tag", "vp", "rho", "name"])
        writer.writerow([10, vp, 1, "medium"])
    (output / "case.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULTS / "fine-fraction-gmsh-results")
    parser.add_argument("--sides", type=float, nargs="+", default=[0, 200, 400, 800, 1200, 1600],
                        help="side of the refined square [m] (0: no patch)")
    parser.add_argument("--h", type=float, nargs="+", default=[40.0], help="element size outside the patch [m]")
    args = parser.parse_args()
    for h in args.h:
        for side in args.sides:
            name = f"h{h:g}_side{side:g}"
            info = generate(args.output / name, side, h)
            print(name, info["elements"], "elements")
