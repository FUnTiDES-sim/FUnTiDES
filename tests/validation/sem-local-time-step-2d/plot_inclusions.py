"""Validate and plot inclusion snapshots, surface gathers, and repeated timings."""
import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/funtides-lts-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


def table(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf8")


def decorate(ax, model, region, receivers=True):
    m = model["mesh"]
    x = (np.arange(m["ex"]) + .5) * m["lx"] / m["ex"]
    z = (np.arange(m["ez"]) + .5) * m["lz"] / m["ez"]
    if region.any():
        ax.contour(x / 1000, z / 1000, region > 0, levels=[.5], colors="black", linewidths=.8)
    if receivers:
        depth = model.get("receivers", {}).get("depth_m", 0)
        label = f"51 récepteurs à {depth:g} m de profondeur" if depth else "51 récepteurs en surface"
        ax.scatter(np.linspace(.5, 5.5, 51), np.full(51, 3-depth/1000), s=8, marker="v", color="#b52e38",
                   clip_on=False, zorder=5, label=label)
    if model.get("boundary_conditions", "").startswith("free_top"):
        ax.plot([0,6],[3,3],color="#286b9c",lw=1,label="Surface libre, pression nulle")
    ax.set(xlim=(0, 6), ylim=(0, 3.08), xlabel="x (km)", ylabel="z (km, vers le haut)")
    ax.set_aspect("equal")


def plot_case(root, gain_power=2., gain_reference_time=1., gain_max=20.):
    model = json.loads((root / "model.json").read_text())
    protocol = json.loads((root / "protocol.json").read_text())
    diagnostic = json.loads((root / "diagnostics.json").read_text())
    spectrum = json.loads((root / "spectral.json").read_text())
    boundary = model.get("boundary_conditions", "rigid")
    absorbing = boundary in ("absorbing", "free_top_absorbing_sides")
    free_top = boundary.startswith("free_top")
    receiver_z = model.get("receivers", {}).get("z", 3000)
    boundary_label = {"absorbing":"Quatre bords absorbants du premier ordre", "rigid":"Bords rigides",
                      "free_top_absorbing_sides":"Surface libre · côtés et fond absorbants",
                      "free_top_rigid_sides":"Surface libre · côtés et fond rigides"}[boundary]
    elements = table(root / "model.csv")
    region = elements["region_id"].reshape(30, 60)
    vp = elements["vp"].reshape(30, 60)
    assert len(elements) == 1800
    assert np.array_equal(elements["element"], np.arange(1800))
    assert all((elements["region_id"] == r["id"]).sum() == r["elements"] for r in model["regions"])
    assert all(np.all(elements["target_dt_s"][elements["region_id"] == r["id"]] == r["target_dt_s"])
               for r in model["regions"])
    assert spectrum["relative_residual"] < .002
    assert diagnostic.get("max_relative_energy_balance_error", diagnostic["max_relative_energy_drift"]) < 1e-4
    if absorbing:
        assert diagnostic["final_energy_fraction"] < 1
    if free_top:
        assert diagnostic["free_surface_max_pressure"] == 0
    snapshots, receivers, history = (table(root / f"{name}.csv") for name in ("snapshots", "receivers", "history"))
    times = np.unique(snapshots["time"])
    assert len(times) == 6 and len(snapshots) == 6 * 7381
    assert len(receivers) == len(history) * 51
    assert np.all(receivers["z"] == receiver_z)
    assert np.all(np.isfinite(snapshots["local"])) and np.all(np.isfinite(snapshots["reference"]))
    assert np.all(np.isfinite(receivers["local"])) and np.all(np.isfinite(receivers["reference"]))
    assert np.allclose(snapshots["local"][snapshots["time"] == 0],
                       snapshots["reference"][snapshots["time"] == 0])
    for t in times:
        assert np.array_equal(snapshots["node"][snapshots["time"] == t], np.arange(7381))
    if free_top:
        surface = snapshots["z"] == 3000
        assert np.all(snapshots["local"][surface] == 0) and np.all(snapshots["reference"][surface] == 0)
    rt = np.unique(receivers["time"])
    xrcv = receivers["x"][:51]
    local = receivers["local"].reshape(-1, 51)
    reference = receivers["reference"].reshape(-1, 51)
    assert np.allclose(receivers["time"].reshape(-1, 51), rt[:, None])
    assert np.allclose(receivers["x"].reshape(-1, 51), xrcv[None, :])
    trace_error = np.linalg.norm(local - reference) / np.linalg.norm(reference)
    assert np.linalg.norm(reference) > 0
    title = {"single_4500": "Petite inclusion à 4 500 m/s dans un fond à 1 500 m/s",
             "single_4500_rigid_surface": "Même inclusion : témoin à bords rigides",
             "single_4500_closed_sides": "Même inclusion : témoin à surface libre et côtés rigides",
             "homogeneous_1500": "Témoin homogène à 1 500 m/s",
             "multiple_noninteger_targets": "Deux inclusions : 4 500 et 3 600 m/s"}[model["name"]]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), layout="constrained")
    image = axes[0].imshow(vp, origin="lower", extent=(0, 6, 0, 3), cmap="viridis", vmin=1500, vmax=4500)
    axes[0].set_title("Modèle de vitesse et acquisition")
    fig.colorbar(image, ax=axes[0], label="Vp (m/s)", shrink=.75)
    fine = model["partition"]["macro_dt_s"] > model["partition"]["cfl"] * 50 / vp
    actual_dt = np.where(fine, model["partition"]["macro_dt_s"] / 2, model["partition"]["macro_dt_s"]) * 1000
    image = axes[1].imshow(actual_dt, origin="lower", extent=(0, 6, 0, 3), cmap="cividis", vmin=6, vmax=12)
    axes[1].set_title(f"Niveau temporel fin : {protocol['fine_elements']} / 1 800 éléments")
    fig.colorbar(image, ax=axes[1], label="Pas du niveau (ms)", shrink=.75)
    for ax in axes:
        decorate(ax, model, region)
        ax.scatter(model["initial_field"]["x"] / 1000, model["initial_field"]["z"] / 1000,
                   marker="*", color="white", edgecolors="black", s=90, label="Centre du champ initial")
    axes[0].legend(loc="lower left", fontsize=8)
    fig.suptitle(title, fontsize=15, fontweight="bold")
    profile = "gaussienne avec image antisymétrique" if free_top else "gaussienne"
    fig.supxlabel(f"{boundary_label} · z = 3 km : surface · Pression initiale {profile}, sans source forcée", fontsize=9)
    fig.savefig(root / "01_modele_et_recepteurs.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(14, 6.5), layout="constrained")
    amplitude = max(np.max(np.abs(snapshots["local"][snapshots["time"] > 0])), .01)
    for ax, time in zip(axes.ravel(), times):
        field = snapshots["local"][snapshots["time"] == time].reshape(61, 121)
        image = ax.imshow(field, origin="lower", extent=(0, 6, 0, 3), cmap="RdBu_r",
                          vmin=-amplitude, vmax=amplitude, interpolation="nearest")
        decorate(ax, model, region)
        ax.set_title(f"t = {time:.2f} s")
    fig.colorbar(image, ax=axes, label="Pression (u.a.), échelle commune", shrink=.8, extend="both")
    fig.suptitle(f"{title}\nLocal 12 / 6 ms · {boundary_label.lower()}", fontsize=15, fontweight="bold")
    fig.supxlabel("Contour noir : inclusion discrétisée · Triangles rouges : récepteurs · Le centre initial à t = 0 est saturé", fontsize=9)
    fig.savefig(root / "02_snapshots.png", dpi=160)
    plt.close(fig)

    difference = local - reference
    # One deterministic time gain, common to every receiver and all three
    # panels. In particular, never normalize the error or each trace separately.
    if not (np.isfinite(gain_power) and gain_power >= 0 and
            np.isfinite(gain_reference_time) and gain_reference_time > 0 and
            np.isfinite(gain_max) and gain_max > 0):
        raise ValueError("Gain power must be nonnegative; reference time and cap must be positive and finite")
    gain = np.minimum((rt / gain_reference_time) ** gain_power, gain_max)
    gained_reference = reference * gain[:, None]
    gained_local = local * gain[:, None]
    gained_difference = difference * gain[:, None]
    amplitude = max(float(np.max(np.abs(gained_reference))), 1e-12)
    common_norm = Normalize(vmin=-amplitude, vmax=amplitude)
    gain_label = f"g(t) = min[(t / {gain_reference_time:g} s)^{gain_power:g}, {gain_max:g}]"
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.6), layout="constrained")
    for ax, values, label in zip(axes, [gained_reference, gained_local, gained_difference],
                                  ["Référence globale : 0,375 ms", "Local : 12 / 6 ms", "Différence local − référence"]):
        image = ax.imshow(values, origin="upper", extent=(.5, 5.5, rt[-1], 0), aspect="auto",
                          cmap="RdBu_r", norm=common_norm, interpolation="nearest")
        ax.set(xlabel="x du récepteur (km)", ylabel="Temps (s)", title=label)
    fig.colorbar(image, ax=axes, label="Pression × gain (u.a.) · échelle commune aux trois panneaux", extend="both")
    fig.suptitle(f"Récepteurs à z = {receiver_z/1000:g} km · {boundary_label.lower()}", fontsize=15, fontweight="bold")
    precision_note = " · Référence : temps float64, forces float32" if protocol.get("reference_double_time_accumulation", False) else ""
    fig.supxlabel(f"{gain_label} · Même gain pour référence, local et erreur ; aucune normalisation par trace\n"
                  f"Erreur L² sur les données sans gain : {trace_error * 100:.3f} %{precision_note}", fontsize=9)
    fig.savefig(root / "03_sismogrammes_surface.png", dpi=160)
    plt.close(fig)
    (root / "display_gain.json").write_text(json.dumps({
        "formula": "min((t/reference_time_s)**power, max_gain)", "power": gain_power,
        "reference_time_s": gain_reference_time, "max_gain": gain_max,
        "actual_max_gain": float(gain.max()), "common_vmin": -amplitude, "common_vmax": amplitude,
        "gain_common_to_reference_local_error": True, "color_scale_common_to_reference_local_error": True,
        "per_trace_normalization": False, "error_metrics_use_raw_data": True,
    }, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(12, 6), layout="constrained")
    for ax, x in zip(axes.ravel(), [1000, 2200, 3000, 4000]):
        receiver = int(np.argmin(abs(xrcv - x)))
        ax.plot(rt, gained_reference[:, receiver], color="#286b9c", label="Global 0,375 ms", lw=1.5)
        ax.plot(rt, gained_local[:, receiver], color="#bf523b", label="Local 12 / 6 ms", ls="--", lw=1.1)
        ax.plot(rt, gained_difference[:, receiver], color="#555555", label="Erreur local − référence", lw=.9)
        ax.set(title=f"x = {xrcv[receiver] / 1000:g} km, z = {receiver_z/1000:g} km", xlabel="Temps (s)", ylabel="Pression × gain (u.a.)")
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=9)
    fig.suptitle(f"Traces avec gain temporel commun et erreur\n{boundary_label}", fontsize=15, fontweight="bold")
    fig.supxlabel(f"{gain_label} · Référence, local et erreur sur les mêmes axes, sans normalisation individuelle", fontsize=9)
    fig.savefig(root / "04_traces_surface.png", dpi=160)
    plt.close(fig)

    rows = table(root / "timing_summary.csv")
    with (root / "timing_samples.csv").open() as stream:
        raw = list(csv.DictReader(stream))
    assert len(rows) == 26 and len(raw) == 26 * protocol["repetitions"]
    for r in rows:
        trials = [s for s in raw if s["case"] == r["case"]]
        assert len(trials) == protocol["repetitions"]
        loop = [float(s["loop_ms"]) for s in trials]
        assert np.isclose(np.median(loop), r["loop_median_ms"])
        assert np.isclose(np.quantile(loop, .25), r["loop_q25_ms"])
        assert np.isclose(np.quantile(loop, .75), r["loop_q75_ms"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
    for is_local, label, color, marker in [(0, "Global", "#286b9c", "o"), (1, "Local", "#078578", "s")]:
        good = rows[(rows["local"] == is_local) & np.isfinite(rows["relative_l2_error"]) & (rows["relative_l2_error"] < .1)]
        axes[0].plot(good["loop_median_ms"], good["relative_l2_error"] * 100, marker + "-", color=color, label=label)
        for r in good:
            axes[0].annotate(f"{r['dt_ms']:.2g} ms", (r["loop_median_ms"], r["relative_l2_error"] * 100),
                             xytext=(4, 4), textcoords="offset points", fontsize=8)
    axes[0].set(xlabel="Boucle médiane (ms)", ylabel="Erreur finale L² (%)", xscale="log", yscale="log",
                title="Coût et précision parmi les pas testés")
    axes[0].legend()
    axes[1].semilogy(history["time"], np.maximum(history["global_macro_amplitude"], 1e-12),
                     color="#b52e38", label="Global 12 ms")
    axes[1].semilogy(history["time"], history["reference_amplitude"], color="#286b9c", label="Référence")
    axes[1].set(xlabel="Temps (s)", ylabel="Amplitude maximale (u.a.)", title="Effet de l’inclusion sur le grand pas global")
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.supxlabel(f"CPU Kokkos {protocol['backend']} · {protocol['repetitions']} répétitions · Préparation et diagnostics exclus", fontsize=9)
    fig.savefig(root / "05_precision_et_stabilite.png", dpi=160)
    plt.close(fig)
    result = {"name": model["name"], "surface_trace_relative_l2_error": float(trace_error),
              "boundary_conditions": model.get("boundary_conditions", "rigid"),
              "final_field_relative_l2_error": float(history["relative_l2_error"][-1]),
              "max_trajectory_relative_l2_error": float(np.max(history["relative_l2_error"]))}
    for budget in [.01, .02, .05]:
        candidates = []
        for is_local in [0, 1]:
            good = rows[(rows["local"] == is_local) & (rows["relative_l2_error"] <= budget)]
            if len(good):
                r = good[np.argmin(good["loop_median_ms"])]
                candidates.append({"local": bool(is_local), "dt_ms": float(r["dt_ms"]),
                                   "loop_ms": float(r["loop_median_ms"]), "error": float(r["relative_l2_error"])})
        result[f"budget_{budget:g}"] = candidates
    (root / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


def compare_surface(root):
    free_top = (root / "single_4500_closed_sides").exists()
    control = "single_4500_closed_sides" if free_top else "single_4500_rigid_surface"
    receiver_z = json.loads((root / "single_4500/model.json").read_text()).get("receivers", {}).get("z",3000)
    folders = [root / name for name in ("single_4500", "homogeneous_1500", control)]
    gathers = [table(p / "receivers.csv") for p in folders]
    histories = [table(p / "history.csv") for p in folders]
    for gather in gathers[1:]:
        assert np.allclose(gather["time"], gathers[0]["time"])
        assert np.array_equal(gather["x"], gathers[0]["x"])
    times = np.unique(gathers[0]["time"])
    fields = [g["reference"].reshape(-1, 51) for g in gathers]
    scattering = fields[0] - fields[1]
    local_scattering = gathers[0]["local"].reshape(-1, 51) - gathers[1]["local"].reshape(-1, 51)
    scattering_error = np.linalg.norm(local_scattering - scattering) / np.linalg.norm(scattering)
    boundary_effect = fields[2] - fields[0]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), layout="constrained")
    for ax, field, title in zip(axes[:2], [scattering, boundary_effect],
                                ["Effet de l’inclusion\nInclusion − fond homogène, bords absorbants",
                                 "Effet des bords\nBords rigides − absorbants, même inclusion"]):
        amplitude = max(np.max(np.abs(field)), 1e-12)
        image = ax.imshow(field, origin="upper", extent=(.5, 5.5, times[-1], 0), aspect="auto", cmap="RdBu_r",
                          vmin=-amplitude, vmax=amplitude, interpolation="nearest")
        ax.set(title=title, xlabel="x du récepteur (km)", ylabel="Temps (s)")
        fig.colorbar(image, ax=ax, label="Différence de pression (u.a.), échelle propre")
    for history, label, color in zip([histories[0], histories[2]], ["Bords absorbants", "Bords rigides"],
                                     ["#078578", "#b52e38"]):
        axes[2].plot(history["time"], history["local_energy"] / history["local_energy"][0], label=label, color=color)
    axes[2].set(xlabel="Temps (s)", ylabel="Énergie discrète / énergie initiale", title="Énergie locale\nMême excitation en surface")
    axes[2].grid(alpha=.2)
    axes[2].legend()
    fig.suptitle("Inclusion, réflexions et absorption" + (" · surface libre" if free_top else " · excitation en surface"), fontsize=15, fontweight="bold")
    fig.supxlabel("Deux panneaux de gauche : références globales à 0,375 ms, temps float64 / forces float32\n"
                  "Énergie : local 12 / 6 ms · Absorption du premier ordre, sans PML", fontsize=9)
    fig.savefig(root / "06_effet_inclusion_et_absorption.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.2), layout="constrained")
    for ax, x in zip(axes.ravel(), [1200, 2200, 3000, 4000]):
        receiver = int(round((x-500)/100))
        ax.plot(times, scattering[:, receiver], color="#286b9c", label="Différence de références globales")
        ax.plot(times, local_scattering[:, receiver], color="#bf523b", ls="--", label="Différence des calculs locaux")
        ax.set(xlim=(1.2, times[-1]), xlabel="Temps (s)", ylabel="Différence de pression (u.a.)",
               title=f"x = {x/1000:g} km, z = {receiver_z/1000:g} km")
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Réponse différentielle de l’inclusion : inclusion − milieu homogène", fontsize=14, fontweight="bold")
    fig.supxlabel(f"Mêmes bords absorbants et même excitation · Amplitudes conservées · Erreur L² différentielle sur toutes les traces : {scattering_error*100:.3f} %", fontsize=9)
    fig.savefig(root / "07_reponse_differentielle_inclusion.png", dpi=160)
    plt.close(fig)
    result = {
        "absorbing_final_energy_fraction": float(histories[0]["local_energy"][-1] / histories[0]["local_energy"][0]),
        "rigid_final_energy_fraction": float(histories[2]["local_energy"][-1] / histories[2]["local_energy"][0]),
        "inclusion_surface_response_l2_norm": float(np.linalg.norm(scattering)),
        "inclusion_surface_response_relative_l2_error": float(scattering_error),
        "rigid_minus_absorbing_surface_difference_l2_norm": float(np.linalg.norm(boundary_effect)),
    }
    assert result["absorbing_final_energy_fraction"] < .5 * result["rigid_final_energy_fraction"]
    assert result["inclusion_surface_response_l2_norm"] > 0
    (root / "surface_comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--compare-surface", action="store_true")
    parser.add_argument("--gain-power", type=float, default=2.)
    parser.add_argument("--gain-reference-time", type=float, default=1., help="Gain reference time in seconds")
    parser.add_argument("--gain-max", type=float, default=20.)
    args = parser.parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    if args.compare_surface:
        compare_surface(args.results)
    else:
        plot_case(args.results, args.gain_power, args.gain_reference_time, args.gain_max)


if __name__ == "__main__":
    main()
