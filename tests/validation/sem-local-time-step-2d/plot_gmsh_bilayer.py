"""Check exported Gmsh benchmark data and plot the actual unstructured mesh.

Run with the numerical Python environment; no Kokkos or GPU import is needed.
"""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/funtides-lts-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np

# Results and reports of the validations
RESULTS = Path(__file__).resolve().parents[3] / "documentation" / "acoustic" / "local-time-step"


def table(root, name):
    return np.atleast_1d(np.genfromtxt(root / f"{name}.csv", delimiter=",", names=True,
                                      dtype=None, encoding="utf8"))


def analyze(root):
    config, protocol, diagnostic, spectrum = (
        json.loads((root / f"{name}.json").read_text())
        for name in ("model", "protocol", "diagnostics", "spectral"))
    nodes, elements, triangles, snapshots, receivers, history, timings, raw = (
        table(root, name) for name in ("nodes", "elements", "triangles", "snapshots",
                                       "receivers", "history", "timing_summary", "timing_samples"))
    assert len(nodes) == protocol["nodes"] and len(elements) == protocol["elements"]
    assert np.array_equal(nodes["node"], np.arange(len(nodes)))
    assert np.array_equal(elements["element"], np.arange(len(elements)))
    assert elements["fine"].sum() == protocol["fine_elements"]
    assert elements["support"].sum() == protocol["support_elements"]
    assert np.all(elements["support"] >= elements["fine"])
    connectivity = np.column_stack([triangles[k] for k in ("n0", "n1", "n2")])
    assert connectivity.min() >= 0 and connectivity.max() < len(nodes)
    tri = mtri.Triangulation(nodes["x"] / 1000, nodes["z"] / 1000, connectivity)
    times = np.unique(snapshots["time"])
    assert len(times) == 6 and len(snapshots) == len(times) * len(nodes)
    top = config["domain"]["lz_m"]
    surface = np.isclose(nodes["z"], top, atol=1e-3, rtol=0)
    for t in times:
        field = snapshots[snapshots["time"] == t]
        assert np.array_equal(field["node"], nodes["node"])
        for k in ("local", "reference"):
            assert np.isfinite(field[k]).all() and np.all(field[k][surface] == 0)
    rt, rx = np.unique(receivers["time"]), np.unique(receivers["receiver"])
    shape = (len(rt), len(rx))
    assert len(receivers) == np.prod(shape) and len(history) == len(rt)
    assert np.allclose(receivers["time"].reshape(shape), rt[:, None])
    assert np.array_equal(receivers["receiver"].reshape(shape), np.broadcast_to(rx, shape))
    xrcv = receivers["x"].reshape(shape)[0]
    assert np.allclose(receivers["x"].reshape(shape), xrcv[None, :])
    assert np.allclose(receivers["z"], top - config["simulation"]["receiver_depth_m"])
    assert np.allclose(history["time"], rt)
    local, reference = (receivers[k].reshape(shape) for k in ("local", "reference"))
    assert np.isfinite(local).all() and np.isfinite(reference).all()
    assert np.linalg.norm(reference) > 0
    assert np.isclose(rt[-1], protocol["duration_s"])
    assert diagnostic["max_relative_energy_balance_error"] < 1e-4
    # Power iteration is an informational estimate, not a certified stability
    # bound. Close dominant eigenvalues can exhaust the iteration budget even
    # when all propagation and energy checks pass. Preserve this status in the
    # analysis instead of silently treating the estimate as converged.
    assert np.isfinite(spectrum["global_dt_limit_estimate_s"]) and spectrum["global_dt_limit_estimate_s"] > 0
    assert np.isfinite(spectrum["relative_residual"]) and spectrum["relative_residual"] >= 0
    e0 = history["energy"][0]
    balance = np.abs(history["energy"] + history["cumulative_dissipation"] - e0) / e0
    assert np.allclose(balance, history["balance_error"], atol=1e-10, rtol=1e-4)
    assert np.isclose(balance.max(), diagnostic["max_relative_energy_balance_error"], atol=1e-10)
    assert np.isclose(history["relative_l2_error"].max(), diagnostic["max_trajectory_relative_l2_error"])
    assert np.all(np.diff(history["cumulative_dissipation"]) >= -1e-12)
    for row in timings:
        samples = raw[raw["case"] == row["case"]]
        assert len(samples) == protocol["repetitions"]
        for field, quantile in (("loop_median_ms", .5), ("loop_q25_ms", .25), ("loop_q75_ms", .75)):
            assert np.isclose(row[field], np.quantile(samples["loop_ms"], quantile))
        assert np.isclose(row["relative_l2_error"], np.median(samples["relative_l2_error"]))
    by_name = {row["case"]: row for row in timings}
    metric = {
        "fine_fraction": float(elements["fine"].mean()),
        "support_fraction": float(elements["support"].mean()),
        "receiver_relative_l2_error_raw": float(np.linalg.norm(local - reference) / np.linalg.norm(reference)),
        "final_energy_fraction": float(history["energy"][-1] / e0),
        "max_relative_energy_balance_error": float(balance.max()),
        "max_trajectory_relative_l2_error": float(history["relative_l2_error"].max()),
        "local_over_global_macro_time": float(by_name["Local"]["loop_median_ms"] /
                                                by_name["Global_macro"]["loop_median_ms"]),
        "local_over_global_fine_time": float(by_name["Local"]["loop_median_ms"] /
                                               by_name["Global_fine"]["loop_median_ms"]),
        "macro_dt_over_global_spectral_estimate": float(protocol["macro_dt_s"] /
                                                        spectrum["global_dt_limit_estimate_s"]),
        "reference_delta_l2": protocol["reference_delta_l2"],
        "spectral_estimate_converged": not spectrum["iteration_limit_reached"],
        "spectral_relative_residual": spectrum["relative_residual"],
        "spectral_iterations": spectrum["iterations"],
        "spectral_note": "Numerical estimate only; convergence status does not certify a stability bound",
        "reference_note": "Same spatial mesh; float64 temporal recurrence, float32 spatial forces",
        "display_gain": "min((time_s / 1.0)**2, 20.0); shared by reference, local and difference",
        "fine_elements_by_vp": {str(v): int(np.sum((elements["vp"] == v) & (elements["fine"] == 1)))
                                 for v in np.unique(elements["vp"])},
    }
    (root / "analysis.json").write_text(json.dumps(metric, indent=2, allow_nan=False) + "\n")

    def decorate(ax):
        inc = config["inclusion"]
        ax.plot(np.array([inc["x_min"], inc["x_max"], inc["x_max"], inc["x_min"], inc["x_min"]]) / 1000,
                np.array([inc["z_min"], inc["z_min"], inc["z_max"], inc["z_max"], inc["z_min"]]) / 1000,
                color="black", lw=.8)
        ax.axhline(config["domain"]["interface_z_m"] / 1000, color="black", lw=.6, ls="--")
        ax.set(xlim=(0, config["domain"]["lx_m"] / 1000), ylim=(0, top / 1000),
               xlabel="x (km)", ylabel="z (km, vers le haut)", aspect="equal")

    vertices = np.stack([np.column_stack([nodes["x"][elements[k]], nodes["z"][elements[k]]])
                         for k in ("n0", "n1", "n2", "n3")], axis=1) / 1000
    # Report an area-based size for the non-square quads; the Gmsh control
    # lengths are targets, not the measured length of every element edge.
    area = .5 * np.abs(np.sum(vertices[:, :, 0] * np.roll(vertices[:, :, 1], -1, axis=1)
                              - vertices[:, :, 1] * np.roll(vertices[:, :, 0], -1, axis=1), axis=1))
    metric["mesh_by_vp"] = {}
    for name in ("upper_layer", "lower_layer", "inclusion"):
        vp = config[name]["vp_m_s"]
        selected = elements["vp"] == vp
        sizes = np.sqrt(area[selected]) * 1000
        metric["mesh_by_vp"][str(vp)] = {
            "elements": int(selected.sum()),
            "target_size_m": config["generation"]["target_final_element_size_m"][name],
            "sqrt_area_median_m": float(np.median(sizes)),
            "sqrt_area_q05_m": float(np.quantile(sizes, .05)),
            "sqrt_area_q95_m": float(np.quantile(sizes, .95)),
        }
    (root / "analysis.json").write_text(json.dumps(metric, indent=2, allow_nan=False) + "\n")
    depth_vertices = vertices.copy()
    depth_vertices[:, :, 1] = top / 1000 - vertices[:, :, 1]
    cmap = ListedColormap(["#e3a84b", "#73b7d0", "#36638d"])
    norm = BoundaryNorm([0, 1000, 2250, 4000], cmap.N)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), layout="constrained", width_ratios=[1.7, 1])
    inc = config["inclusion"]
    for ax in axes:
        collection = PolyCollection(depth_vertices, array=elements["vp"], cmap=cmap, norm=norm,
                                    edgecolors="#263746", linewidths=.35)
        ax.add_collection(collection)
        ax.set(xlabel="x (km)", ylabel="Profondeur (km)", aspect="equal")
        ax.axhline((top - config["domain"]["interface_z_m"]) / 1000, color="black", lw=1)
        ax.plot(np.array([inc["x_min"], inc["x_max"], inc["x_max"], inc["x_min"], inc["x_min"]]) / 1000,
                (top - np.array([inc["z_min"], inc["z_min"], inc["z_max"], inc["z_max"], inc["z_min"]])) / 1000,
                color="black", lw=1.2)
    axes[0].set(xlim=(0, config["domain"]["lx_m"] / 1000), ylim=(top / 1000, 0),
                title=f"Modèle complet · {len(elements):,} quadrilatères SEM d’ordre 2".replace(",", " "))
    for name, depth in (("upper_layer", .65), ("lower_layer", 2.65)):
        v = config[name]["vp_m_s"]
        target = metric["mesh_by_vp"][str(v)]["target_size_m"]
        axes[0].text(.18, depth, f"{v:,} m/s · cible {target:g} m".replace(",", " "),
                     bbox=dict(facecolor="white", alpha=.9, edgecolor="none"), fontsize=10)
    axes[1].set(xlim=(inc["x_min"] / 1000 - .35, inc["x_max"] / 1000 + .35),
                ylim=((top - inc["z_min"]) / 1000 + .3, (top - inc["z_max"]) / 1000 - .3),
                title="Zoom : raffinement autour de l’inclusion")
    inclusion_stats = metric["mesh_by_vp"][str(config["inclusion"]["vp_m_s"])]
    axes[1].annotate(f"{config['inclusion']['vp_m_s']:g} m/s · cible {inclusion_stats['target_size_m']:g} m\n"
                     f"{inclusion_stats['elements']} éléments",
                     xy=((inc["x_min"] + inc["x_max"]) / 2000, (top - inc["z_min"]) / 1000 - .04),
                     xytext=((inc["x_min"] + inc["x_max"]) / 2000, (top - inc["z_min"]) / 1000 + .18),
                     ha="center", va="center", fontsize=10, arrowprops=dict(arrowstyle="->", color="black"),
                     bbox=dict(facecolor="white", alpha=.95, edgecolor="none"))
    fig.colorbar(collection, ax=axes, ticks=[600, 1500, 3000], label="Vitesse Vp (m/s)", shrink=.75)
    fig.suptitle("Bicouche Gmsh avec inclusion lente à 600 m/s", fontsize=16)
    fig.supxlabel("Traits : vrais éléments du fichier Gmsh · tailles cibles : 100 m / 200 m / 40 m\n"
                  "L’inclusion de 400 × 400 m se trouve entre 1,9 et 2,3 km de profondeur.", fontsize=10)
    fig.savefig(root / "00_modele_bicouche_gmsh.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(11, 9), layout="constrained")
    for ax, values, title, label, cmap in (
        (axes[0], elements["vp"], f"Maillage adapté à la vitesse : {len(elements):,} quadrilatères".replace(",", " "), "Vp (m/s)", "viridis"),
        (axes[1], np.where(elements["fine"], protocol["fine_dt_s"], protocol["macro_dt_s"]) * 1000,
         f"Partition temporelle : {metric['fine_fraction']:.2%} fins, {metric['support_fraction']:.2%} de support",
         "Pas du niveau (ms)", "cividis")):
        collection = PolyCollection(vertices, array=values, cmap=cmap, edgecolors="#555555", linewidths=.18)
        ax.add_collection(collection)
        fig.colorbar(collection, ax=ax, label=label, shrink=.8)
        decorate(ax)
        ax.scatter(xrcv / 1000, np.full(len(rx), (top - config["simulation"]["receiver_depth_m"]) / 1000),
                   marker="v", s=9, color="#b52e38", label="Récepteurs à 50 m")
        ax.scatter(config["simulation"]["source_x_m"] / 1000,
                   (top - config["simulation"]["source_depth_m"]) / 1000,
                   marker="*", s=110, c="white", edgecolors="black", label="Excitation initiale")
        ax.set_title(title)
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Bicouche Gmsh · surface libre, côtés et fond absorbants")
    fig.savefig(root / "01_maillage_et_partition.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), layout="constrained")
    amplitude = np.max(np.abs(snapshots["local"][snapshots["time"] > 0]))
    for ax, time in zip(axes.flat, times):
        field = snapshots["local"][snapshots["time"] == time]
        image = ax.tripcolor(tri, field, shading="gouraud", cmap="RdBu_r", vmin=-amplitude, vmax=amplitude)
        decorate(ax)
        ax.set_title(f"t = {time:.2f} s")
    fig.colorbar(image, ax=axes, label="Pression (u.a.) · échelle commune", shrink=.8, extend="both")
    fig.suptitle("Propagation locale sur le maillage Gmsh · pic initial saturé")
    fig.savefig(root / "02_propagation.png", dpi=150)
    plt.close(fig)

    gain = np.minimum(rt ** 2, 20)[:, None]
    difference = local - reference
    amplitude = np.max(np.abs(reference * gain))
    fig, axes = plt.subplots(1, 3, figsize=(14, 6), layout="constrained")
    for ax, values, label in zip(axes, (reference, local, difference), ("Référence", "Local", "Local − référence")):
        image = ax.imshow(values * gain, origin="upper", aspect="auto", interpolation="nearest",
                          extent=(xrcv[0] / 1000, xrcv[-1] / 1000, rt[-1], rt[0]),
                          cmap="RdBu_r", vmin=-amplitude, vmax=amplitude)
        ax.set(xlabel="x du récepteur (km)", ylabel="Temps (s)", title=label)
    fig.colorbar(image, ax=axes, label="Pression × gain (u.a.) · même échelle")
    fig.suptitle(f"Récepteurs à 50 m · erreur L² sans gain : {metric['receiver_relative_l2_error_raw']:.3%}")
    fig.supxlabel("Gain commun g(t) = min[(t / 1 s)², 20] · aucune normalisation par trace", fontsize=10)
    fig.savefig(root / "03_sismogrammes.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), layout="constrained")
    for ax, x in zip(axes.flat, (1000, 2200, 3000, 4000)):
        j = np.argmin(abs(xrcv - x))
        for values, label, color, style in ((reference, "Référence", "#286b9c", "-"),
                                            (local, "Local", "#bf523b", "--"),
                                            (difference, "Erreur", "#555555", "-")):
            ax.plot(rt, values[:, j] * gain[:, 0], label=label, color=color, ls=style, lw=1)
        ax.set(xlabel="Temps (s)", ylabel="Pression × gain (u.a.)", title=f"x = {xrcv[j] / 1000:g} km")
        ax.grid(alpha=.2)
    axes[0, 0].legend()
    fig.suptitle("Traces et erreurs · même gain temporel et mêmes axes")
    fig.savefig(root / "04_traces.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    labels = ["Global petit pas", "Global grand pas", "Local"]
    med = timings["loop_median_ms"]
    axes[0, 0].bar(labels, med, color=["#286b9c", "#6999b9", "#078578"])
    axes[0, 0].errorbar(np.arange(3), med,
                       yerr=[med - timings["loop_q25_ms"], timings["loop_q75_ms"] - med], fmt="none", c="black")
    axes[0, 0].set(ylabel="Boucle médiane et IQR (ms)",
                   title=f"Temps CPU enregistrés · {protocol['repetitions']} répétitions")
    axes[0, 1].bar(labels, timings["relative_l2_error"] * 100, color=["#286b9c", "#6999b9", "#078578"])
    axes[0, 1].set(ylabel="Erreur finale L² pondérée par la masse (%)", title="Même référence spatiale et temporelle")
    axes[1, 0].plot(rt, history["energy"] / e0, label="E / E₀")
    axes[1, 0].plot(rt, history["cumulative_dissipation"] / e0, label="Dissipation cumulée / E₀")
    axes[1, 0].set(xlabel="Temps (s)", title="Bilan d’énergie local")
    axes[1, 0].legend()
    axes[1, 1].semilogy(rt, np.maximum(balance, 1e-16), label="Erreur de bilan énergétique")
    axes[1, 1].semilogy(rt, np.maximum(history["relative_l2_error"], 1e-16), label="Erreur du champ vs référence")
    axes[1, 1].set(xlabel="Temps (s)", ylabel="Erreur relative", title="Contrôles sur toute la trajectoire")
    axes[1, 1].legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("Bicouche Gmsh · précision, coût et dissipation")
    fig.savefig(root / "05_precision_et_cout.png", dpi=150)
    plt.close(fig)
    print(json.dumps(metric, indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS / "gmsh-bilayer-results")
    analyze(parser.parse_args().results)
