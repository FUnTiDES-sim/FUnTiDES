"""Convert a linear Q4/Hex8 Gmsh file to the native reader's MSH 2.2 ASCII format."""
import argparse
from pathlib import Path
import gmsh


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path);p.add_argument("output",type=Path)
    p.add_argument("--plane",choices=["xz","xy"],default="xz")
    args=p.parse_args()
    if args.input.resolve()==args.output.resolve():raise ValueError("Use a distinct output path")
    gmsh.initialize(["gmsh","-nopopup"])
    try:
        gmsh.open(str(args.input))
        dim = gmsh.model.getDimension()
        expected = 5 if dim == 3 else 3
        if dim not in (2, 3) or list(gmsh.model.mesh.getElementTypes(dim)) != [expected]:
            raise ValueError("Expected only linear Q4 quadrilaterals (2D) or Hex8 hexahedra (3D); triangles and curved/high-order geometry are not silently converted")
        if dim == 3 and args.plane == "xy":
            raise ValueError("--plane xy applies only to 2D meshes")
        if dim == 2 and args.plane=="xy":
            gmsh.model.mesh.affineTransform([1,0,0,0, 0,0,-1,0, 0,1,0,0, 0,0,0,1])
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
        gmsh.option.setNumber("Mesh.Binary",0)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        gmsh.write(str(args.output))
    finally:gmsh.finalize()


if __name__=="__main__":main()
