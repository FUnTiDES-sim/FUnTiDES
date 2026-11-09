"""Inspect and display a Gmsh mesh without running any wave propagation."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/funtides-lts-matplotlib")
import gmsh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np


def outside_cells(root):
    """Canonical exterior cell geometry, independent of Gmsh node numbering."""
    config = json.loads((root / "model.json").read_text())
    gmsh.initialize(["gmsh", "-nopopup"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(root / "bilayer.msh"))
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        coords = {int(t): tuple(np.round(p[[0, 2]], 6)) for t, p in zip(tags, np.asarray(xyz).reshape(-1, 3))}
        result = []
        for name in ("upper_layer", "lower_layer"):
            for entity in gmsh.model.getEntitiesForPhysicalGroup(2, config[name]["physical_tag"]):
                _, _, connectivity = gmsh.model.mesh.getElements(2, int(entity))
                result.extend((name, tuple(sorted(coords[int(t)] for t in cell)))
                              for cell in np.asarray(connectivity[0]).reshape(-1, 4))
        return sorted(result)
    finally:
        gmsh.finalize()


def preview(root):
    config = json.loads((root / "model.json").read_text())
    gmsh.initialize(["gmsh", "-nopopup"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(root / "bilayer.msh"))
        tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.asarray(coords).reshape(-1, 3)
        index = {int(tag): i for i, tag in enumerate(tags)}
        cells, speeds, regions = [], [], []
        for name in ("upper_layer", "lower_layer", "inclusion"):
            material = config[name]
            for entity in gmsh.model.getEntitiesForPhysicalGroup(2, material["physical_tag"]):
                types, element_tags, node_tags = gmsh.model.mesh.getElements(2, int(entity))
                if list(types) != [3]:
                    raise ValueError("Preview expects only four-node quadrilaterals")
                connectivity = np.asarray(node_tags[0]).reshape(-1, 4)
                cells.extend([[index[int(tag)] for tag in row] for row in connectivity])
                speeds.extend([material["vp_m_s"]] * len(connectivity))
                regions.extend([name] * len(connectivity))
        cells, speeds, regions = np.asarray(cells), np.asarray(speeds), np.asarray(regions)
    finally:
        gmsh.finalize()
    vertices = coords[cells][:, :, [0, 2]]
    edges = np.roll(vertices, -1, axis=1) - vertices
    rectangles = np.all((np.abs(edges[:, :, 0]) < 1e-6) | (np.abs(edges[:, :, 1]) < 1e-6), axis=1)
    areas = .5 * np.abs(np.sum(vertices[:, :, 0] * np.roll(vertices[:, :, 1], -1, axis=1)
                               - vertices[:, :, 1] * np.roll(vertices[:, :, 0], -1, axis=1), axis=1))
    domain = config["domain"]
    assert np.all(areas > 0) and np.isclose(areas.sum(), domain["lx_m"] * domain["lz_m"])
    edge_counts = Counter(tuple(sorted((int(a), int(b)))) for cell in cells
                          for a, b in zip(cell, np.roll(cell, -1)))
    assert all(n in (1, 2) for n in edge_counts.values())
    for (a, b), count in edge_counts.items():
        if count == 1:
            p, q = coords[a], coords[b]
            assert (np.isclose(p[0], q[0]) and (np.isclose(p[0], 0) or np.isclose(p[0], domain["lx_m"]))) or (
                np.isclose(p[2], q[2]) and (np.isclose(p[2], 0) or np.isclose(p[2], domain["lz_m"])))
    assert rectangles[regions == "upper_layer"].all()
    centers = vertices.mean(axis=1)
    refinement = config["mesh"].get("interior_refinement")
    inclusion_regular = regions == "inclusion"
    if refinement:
        inclusion_regular &= ((centers[:, 0] > refinement["x_min"]) & (centers[:, 0] < refinement["x_max"]) &
                              (centers[:, 1] > refinement["z_min"]) & (centers[:, 1] < refinement["z_max"]))
    assert rectangles[inclusion_regular].all()
    patch = config["mesh"]["inclusion_patch"]
    ztransition = config["mesh"]["interface_transition_bottom_z_m"]
    lower_regular = (regions == "lower_layer") & (centers[:, 1] < ztransition) & (
        (centers[:, 0] < patch["x_min"]) | (centers[:, 0] > patch["x_max"]) |
        (centers[:, 1] < patch["z_min"]))
    assert rectangles[lower_regular].all()
    sizes = np.ptp(vertices, axis=1)
    checks = dict(elements=len(cells), vertex_nodes=len(coords), rectangular_elements=int(rectangles.sum()),
                  boundary_edges=sum(n == 1 for n in edge_counts.values()),
                  interior_edges=sum(n == 2 for n in edge_counts.values()),
                  area_m2=float(areas.sum()), no_unmatched_interior_edges=True, propagation_run=False,
                  regions={})
    if refinement:
        assert outside_cells(root) == outside_cells(root.parent / "gmsh-structured-results")
        checks["outside_mesh_identical_to_baseline_at_1e_minus6_m"] = True
    for name in ("upper_layer", "lower_layer", "inclusion"):
        mask = regions == name
        assert mask.sum() == config["generation"]["element_counts"][name]
        regular = lower_regular if name == "lower_layer" else inclusion_regular if name == "inclusion" else mask
        checks["regions"][name] = dict(elements=int(mask.sum()), regular_block_elements=int(regular.sum()),
                                       regular_dx_m=float(np.median(sizes[regular, 0])),
                                       regular_dz_m=float(np.median(sizes[regular, 1])))
    (root / "mesh_checks.json").write_text(json.dumps(checks, indent=2) + "\n")

    depth_vertices = vertices / 1000
    depth_vertices[:, :, 1] = domain["lz_m"] / 1000 - depth_vertices[:, :, 1]
    cmap = ListedColormap(["#e3a84b", "#73b7d0", "#36638d"])
    norm = BoundaryNorm([0, 1000, 2250, 4000], cmap.N)
    fig = plt.figure(figsize=(15, 8), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.6, 1])
    full, interface_ax, inclusion_ax = fig.add_subplot(grid[:, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 1])
    axes = [full, interface_ax, inclusion_ax]
    interface_depth = (domain["lz_m"] - domain["interface_z_m"]) / 1000
    transition_depth = (domain["lz_m"] - ztransition) / 1000
    inc = config["inclusion"]
    for ax in axes:
        collection = PolyCollection(depth_vertices, array=speeds, cmap=cmap, norm=norm,
                                    edgecolors="#213341", linewidths=.55)
        ax.add_collection(collection)
        ax.axhline(interface_depth, color="black", lw=1.4)
        ax.axhline(transition_depth, color="#555555", lw=.9, ls="--")
        ax.set(xlabel="x (km)", ylabel="Profondeur (km)", aspect="equal")
    full.set(xlim=(0, domain["lx_m"] / 1000), ylim=(domain["lz_m"] / 1000, 0),
             title=f"Maillage complet · {len(cells):,} quadrilatères conformes".replace(",", " "))
    for depth, text in ((.55, "Couche 1 : 1 500 m/s\nRectangles 100 × 93,75 m"),
                        (2.6, "Couche 2 : 3 000 m/s\nCarrés 200 × 200 m")):
        full.text(.15, depth, text, fontsize=10, bbox=dict(facecolor="white", alpha=.95, edgecolor="none"))
    interface_ax.set(xlim=(1.8, 4.2), ylim=(1.95, 1.25), title="Raccord 100 → 200 m sous l’interface")
    if refinement:
        inclusion_ax.set(xlim=(2.94, 3.06), ylim=(2.16, 2.04),
                         title=f"Intérieur de l’inclusion : carrés de {refinement['element_size_m']:g} m")
    else:
        inclusion_ax.set(xlim=(2.3, 3.7), ylim=(2.7, 1.8), title="Inclusion à 600 m/s : grille régulière 40 × 40 m")
        inclusion_ax.annotate("100 carrés dans l’inclusion", xy=(3, 2.26), xytext=(3, 2.58), ha="center", fontsize=10,
                              arrowprops=dict(arrowstyle="->"), bbox=dict(facecolor="white", alpha=.95, edgecolor="none"))
    fig.colorbar(collection, ax=axes, ticks=[600, 1500, 3000], label="Vp (m/s)", shrink=.65)
    fig.suptitle("Bicouche Gmsh · blocs rectangulaires et raccords locaux", fontsize=17)
    fig.supxlabel("Maillage seul, avant modélisation · raffinement intérieur à 600 m/s, extérieur inchangé" if refinement else
                  "Maillage seul, avant modélisation · irrégularités limitées au raccord des couches et au pourtour de l’inclusion",
                  fontsize=11)
    fig.savefig(root / "00_maillage_structure.png", dpi=180)
    plt.close(fig)
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path(__file__).parent / "gmsh-structured-results")
    preview(parser.parse_args().results)
