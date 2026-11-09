"""Generate a conforming Hex8 box for the FUnTiDES acoustic SEM importer."""
import argparse
from pathlib import Path

import gmsh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cells", type=int, default=2, help="Cells along each axis")
    parser.add_argument("--msh-version", type=float, choices=[2.2, 4.1], default=2.2)
    args = parser.parse_args()
    if args.cells < 1:
        parser.error("--cells must be positive")
    gmsh.initialize(["gmsh", "-nopopup"])
    try:
        gmsh.model.add("hex_box")
        geo = gmsh.model.geo
        points = [geo.addPoint(x, y, 0) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
        edges = [geo.addLine(points[i], points[(i + 1) % 4]) for i in range(4)]
        surface = geo.addPlaneSurface([geo.addCurveLoop(edges)])
        for edge in edges:
            geo.mesh.setTransfiniteCurve(edge, args.cells + 1)
        geo.mesh.setTransfiniteSurface(surface)
        geo.mesh.setRecombine(2, surface)
        extrusion = geo.extrude([(2, surface)], 0, 0, 1, [args.cells], recombine=True)
        geo.synchronize()
        volumes = [tag for dim, tag in extrusion if dim == 3]
        top = extrusion[0][1]
        sides = [tag for dim, tag in extrusion[2:] if dim == 2]
        gmsh.model.addPhysicalGroup(3, volumes, 10)
        gmsh.model.addPhysicalGroup(2, [top], 1)
        gmsh.model.addPhysicalGroup(2, [surface] + sides, 2)
        gmsh.model.mesh.generate(3)
        if list(gmsh.model.mesh.getElementTypes(3)) != [5]:
            raise RuntimeError("Expected only Hex8 cells")
        gmsh.option.setNumber("Mesh.MshFileVersion", args.msh_version)
        gmsh.option.setNumber("Mesh.Binary", 0)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(args.output))
        args.output.with_suffix(".csv").write_text("physical_tag,vp,rho\n10,1500,1000\n")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
