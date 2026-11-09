"""Legacy all-unstructured Gmsh mesh, kept to reproduce the initial benchmark."""
import argparse
import csv
import json
from pathlib import Path
import gmsh
import numpy as np

# Results and reports of the validations
RESULTS = Path(__file__).resolve().parents[3] / "documentation" / "acoustic" / "local-time-step"


def main():
    here = Path(__file__).parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=here / "gmsh_bilayer.json")
    parser.add_argument("--output", type=Path, default=RESULTS / "gmsh-bilayer-results")
    args = parser.parse_args()
    c = json.loads(args.config.read_text())
    domain, inc, mesh = c["domain"], c["inclusion"], c["mesh"]
    lx, lz, interface = (float(domain[k]) for k in ("lx_m", "lz_m", "interface_z_m"))
    assert 0 < interface < lz
    assert 0 < inc["x_min"] < inc["x_max"] < lx
    assert 0 < inc["z_min"] < inc["z_max"] < interface
    assert mesh["design_fmax_hz"] > 0 and mesh["elements_per_wavelength"] > 0
    args.output.mkdir(parents=True, exist_ok=True)
    materials = [("upper_layer", c["upper_layer"]), ("lower_layer", c["lower_layer"]), ("inclusion", inc)]
    targets = {name: mat["vp_m_s"] / (mesh["design_fmax_hz"] * mesh["elements_per_wavelength"])
               for name, mat in materials}
    assert all(mat["vp_m_s"] > 0 and mat["rho"] > 0 for _, mat in materials)
    assert len({mat["physical_tag"] for _, mat in materials}) == 3
    gmsh.initialize(["gmsh", "-nopopup"])
    try:
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.model.add("velocity_adapted_bilayer")
        geo = gmsh.model.geo
        coords = [(0,0),(lx,0),(lx,interface),(lx,lz),(0,lz),(0,interface),
                  (inc["x_min"],inc["z_min"]),(inc["x_max"],inc["z_min"]),
                  (inc["x_max"],inc["z_max"]),(inc["x_min"],inc["z_max"])]
        points = [geo.addPoint(x,0,z) for x,z in coords]
        exterior = [geo.addLine(points[i],points[(i+1)%6]) for i in range(6)]
        interface_line = geo.addLine(points[5],points[2])
        inner = [geo.addLine(points[6+i],points[6+(i+1)%4]) for i in range(4)]
        upper_loop = geo.addCurveLoop([interface_line,exterior[2],exterior[3],exterior[4]])
        lower_loop = geo.addCurveLoop([exterior[0],exterior[1],-interface_line,exterior[5]])
        hole_loop = geo.addCurveLoop([-e for e in reversed(inner)])
        inclusion_loop = geo.addCurveLoop(inner)
        surfaces = [geo.addPlaneSurface([upper_loop]),geo.addPlaneSurface([lower_loop,hole_loop]),
                    geo.addPlaneSurface([inclusion_loop])]
        geo.synchronize()
        for surface, (name,mat) in zip(surfaces,materials):
            gmsh.model.addPhysicalGroup(2,[surface],int(mat["physical_tag"]))
            gmsh.model.setPhysicalName(2,int(mat["physical_tag"]),name)
        gmsh.model.addPhysicalGroup(1,[exterior[3]],1)
        gmsh.model.setPhysicalName(1,1,"free_surface")
        gmsh.model.addPhysicalGroup(1,[e for i,e in enumerate(exterior) if i!=3],2)
        gmsh.model.setPhysicalName(1,2,"absorbing_sides_bottom")
        # All-triangle subdivision produces conforming quads without hanging
        # nodes. The control length is doubled because subdivision adds midpoints.
        def box(x0,x1,z0,z1,value,outside):
            f=gmsh.model.mesh.field.add("Box")
            for key,v in {"XMin":x0,"XMax":x1,"YMin":-1,"YMax":1,"ZMin":z0,"ZMax":z1,
                          "VIn":value,"VOut":outside,"Thickness":mesh["transition_m"]}.items():
                gmsh.model.mesh.field.setNumber(f,key,v)
            return f
        coarse=2*targets["lower_layer"]
        upper=box(-lx,2*lx,interface,2*lz,2*targets["upper_layer"],coarse)
        inclusion=box(inc["x_min"],inc["x_max"],inc["z_min"],inc["z_max"],2*targets["inclusion"],coarse)
        f=gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(f,"FieldsList",[upper,inclusion])
        gmsh.model.mesh.field.setAsBackgroundMesh(f)
        for key,value in {"Mesh.MeshSizeFromPoints":0,"Mesh.MeshSizeFromCurvature":0,
                          "Mesh.MeshSizeExtendFromBoundary":0,"Mesh.Algorithm":6,
                          "Mesh.SubdivisionAlgorithm":1,"Mesh.ElementOrder":1,
                          "Mesh.MshFileVersion":2.2,"Mesh.Binary":0,"Mesh.SaveAll":0}.items():
            gmsh.option.setNumber(key,value)
        gmsh.model.mesh.generate(2)
        types,tags,_ = gmsh.model.mesh.getElements(2)
        assert list(types)==[3], f"Expected only 4-node quadrilaterals, got {types}"
        quality=np.asarray(gmsh.model.mesh.getElementQualities(tags[0]))
        assert np.isfinite(quality).all() and quality.min()>0
        counts={}
        for surface,(name,mat) in zip(surfaces,materials):
            _, element_tags, _ = gmsh.model.mesh.getElements(2,surface)
            counts[name]=sum(len(t) for t in element_tags)
        gmsh.write(str(args.output / "bilayer.msh"))
        gmsh.write(str(args.output / "bilayer.geo_unrolled"))
        c["generation"]={"gmsh_version":gmsh.__version__,"format":"MSH 2.2 ASCII, Q4, XZ plane",
                         "target_final_element_size_m":targets,"element_counts":counts,
                         "minimum_quality_minSICN":float(quality.min()),
                         "method":"velocity-based size field; conforming triangle-to-quad subdivision"}
    finally:
        gmsh.finalize()
    (args.output / "model.json").write_text(json.dumps(c,indent=2)+"\n")
    with (args.output / "materials.csv").open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["physical_tag","vp","rho","name"])
        for name,mat in materials:writer.writerow([mat["physical_tag"],mat["vp_m_s"],mat["rho"],name])
    sim=c["simulation"]
    (args.output / "simulation.txt").write_text(
        f"{sim['duration_s']} {sim['source_x_m']} {lz-sim['source_depth_m']} {sim['source_sigma_m']} "
        f"{sim['receiver_depth_m']} {sim['cfl']}\n")
    print(json.dumps(c["generation"],indent=2))


if __name__=="__main__":
    main()
