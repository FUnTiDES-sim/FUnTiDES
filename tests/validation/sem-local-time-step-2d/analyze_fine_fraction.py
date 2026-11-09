"""Validate the fraction benchmark CSVs and plot the measured accuracy/cost frontier."""

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/funtides-lts-matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Results and reports of the validations
RESULTS = Path(__file__).resolve().parents[3] / "documentation" / "acoustic" / "local-time-step"


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS / "fraction-results")
    args = parser.parse_args()
    root = args.results
    protocol = json.loads((root / "fraction_protocol.json").read_text())
    scenarios = read_csv(root / "fraction_scenarios.csv")
    summary = read_csv(root / "fraction_summary.csv")
    samples = read_csv(root / "fraction_samples.csv")
    assert len(scenarios) == 8
    assert len(summary) == 8 * 12
    assert len(samples) == len(summary) * protocol["repetitions"]
    keys = {(r["scenario"], r["case"]) for r in summary}
    assert len(keys) == len(summary)
    assert {(r["scenario"], r["case"]) for r in samples} == keys
    for r in summary:
        trials = [s for s in samples if (s["scenario"], s["case"]) == (r["scenario"], r["case"])]
        assert {int(s["trial"]) for s in trials} == set(range(protocol["repetitions"]))
        assert np.isclose(float(r["dt_ms"]) * int(r["steps"]), 1800)
        setup = np.array([float(s["setup_ms"]) for s in trials])
        loop = np.array([float(s["loop_ms"]) for s in trials])
        assert np.all(setup > 0) and np.all(loop > 0)
        expected = {
            "setup_median_ms": np.median(setup), "loop_median_ms": np.median(loop),
            "loop_q25_ms": np.quantile(loop, .25), "loop_q75_ms": np.quantile(loop, .75),
            "total_median_ms": np.median(setup + loop),
        }
        for key, value in expected.items():
            assert np.isclose(float(r[key]), value, rtol=1e-9), (r["case"], key)
        errors = np.array([float(s["relative_l2_error"]) for s in trials])
        assert all(bool(int(s["finite"])) == np.isfinite(e) for s, e in zip(trials, errors))
        assert np.isclose(float(r["relative_l2_error"]), np.median(errors))
        assert bool(int(r["finite"])) == np.isfinite(errors).all()
    for r in scenarios:
        if r["scenario"].startswith("fast_rows_"):
            rows = int(r["scenario"].split("_")[-1])
            assert int(r["fine_elements"]) == rows * 60
        # Reference refinement discrepancy must stay below 0.1%, well below the
        # smallest 0.5% error budget used for the performance comparisons.
        assert float(r["reference_delta_l2"]) < 1e-3, r

    def valid(rows):
        return [r for r in rows if r["finite"] == "1" and np.isfinite(float(r["relative_l2_error"]))]

    def best(rows, local, budget):
        eligible = [r for r in valid(rows) if bool(int(r["local"])) == local
                    and float(r["relative_l2_error"]) <= budget]
        return min(eligible, key=lambda r: float(r["loop_median_ms"]), default=None)

    frontier = []
    for scenario in scenarios:
        name = scenario["scenario"]
        group = [r for r in summary if r["scenario"] == name]
        for budget in [.005, .01, .02, .05]:
            local, glob = best(group, True, budget), best(group, False, budget)
            if local is None or glob is None:
                continue
            frontier.append({
                "scenario": name, "fine_fraction": int(scenario["fine_elements"]) / int(scenario["elements"]),
                "budget_percent": budget * 100, "local_case": local["case"], "global_case": glob["case"],
                "local_error_percent": float(local["relative_l2_error"]) * 100,
                "global_error_percent": float(glob["relative_l2_error"]) * 100,
                "local_ms": float(local["loop_median_ms"]), "global_ms": float(glob["loop_median_ms"]),
                "speedup_global_over_local": float(glob["loop_median_ms"]) / float(local["loop_median_ms"]),
            })
    with (root / "fraction_frontier.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(frontier[0]))
        writer.writeheader()
        writer.writerows(frontier)

    refined = []
    if (root / "budget_protocol.json").exists():
        budget_protocol = json.loads((root / "budget_protocol.json").read_text())
        budget_rows = read_csv(root / "budget_summary.csv")
        budget_samples = read_csv(root / "budget_samples.csv")
        search = read_csv(root / "budget_search.csv")
        assert len(budget_rows) == 16
        assert len(budget_samples) == 16 * budget_protocol["repetitions"]
        assert {(r["scenario"], r["local"]) for r in budget_rows} == {
            (s["scenario"], str(local)) for s in scenarios for local in (0, 1)}
        for r in budget_rows:
            key = r["scenario"], r["local"]
            count = int(r["steps"])
            assert np.isclose(float(r["dt_ms"]) * count, 1800)
            assert float(r["relative_l2_error"]) <= budget_protocol["budget"]
            evaluations = {int(s["steps"]): float(s["relative_l2_error"]) for s in search
                           if (s["scenario"], s["local"]) == key}
            assert evaluations[count] <= budget_protocol["budget"]
            if count != budget_protocol["min_steps"]:
                assert evaluations[count - 1] > budget_protocol["budget"]
            trials = [s for s in budget_samples if (s["scenario"], s["local"]) == key]
            assert {int(s["trial"]) for s in trials} == set(range(budget_protocol["repetitions"]))
            assert all(int(s["steps"]) == count and float(s["relative_l2_error"]) <= budget_protocol["budget"]
                       for s in trials)
            setup = np.array([float(s["setup_ms"]) for s in trials])
            loop = np.array([float(s["loop_ms"]) for s in trials])
            assert np.all(setup > 0) and np.all(loop > 0)
            for column, value in {
                "setup_median_ms": np.median(setup), "loop_median_ms": np.median(loop),
                "loop_q25_ms": np.quantile(loop, .25), "loop_q75_ms": np.quantile(loop, .75),
                "total_median_ms": np.median(setup + loop),
                "relative_l2_error": np.median([float(s["relative_l2_error"]) for s in trials]),
            }.items():
                assert np.isclose(float(r[column]), value, rtol=1e-9)
        for s in scenarios:
            group = [r for r in budget_rows if r["scenario"] == s["scenario"]]
            local = next(r for r in group if r["local"] == "1")
            glob = next(r for r in group if r["local"] == "0")
            refined.append({
                "scenario": s["scenario"], "fine_fraction": int(s["fine_elements"]) / int(s["elements"]),
                "local_dt_ms": float(local["dt_ms"]), "global_dt_ms": float(glob["dt_ms"]),
                "local_error_percent": float(local["relative_l2_error"]) * 100,
                "global_error_percent": float(glob["relative_l2_error"]) * 100,
                "local_ms": float(local["loop_median_ms"]), "global_ms": float(glob["loop_median_ms"]),
                "speedup_global_over_local": float(glob["loop_median_ms"]) / float(local["loop_median_ms"]),
            })
        with (root / "budget_comparison.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(refined[0]))
            writer.writeheader()
            writer.writerows(refined)
        print(f"Validated {len(budget_samples)} refinement samples and {len(search)} search evaluations.")
        for r in refined:
            print(f"Refined {r['scenario']}: local {r['local_dt_ms']:.4f} ms, global {r['global_dt_ms']:.4f} ms; "
                  f"global/local = {r['speedup_global_over_local']:.3f}")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.3), layout="constrained")
    fig.suptitle("SEM 2D : proportion au petit pas, coût et précision", fontsize=16, fontweight="bold")
    for ax, name, title in zip(axes[:2], ["fast_rows_3", "original_three_layers"],
                               ["Couche rapide : 10 % des éléments", "Trois couches initiales : 66,7 % au petit pas"]):
        rows = valid([r for r in summary if r["scenario"] == name])
        for local, label, color, marker in [(False, "Global", "#286b9c", "o"), (True, "Local", "#078578", "s")]:
            points = [r for r in rows if bool(int(r["local"])) == local and float(r["relative_l2_error"]) < .2]
            x = np.array([float(r["loop_median_ms"]) for r in points])
            y = np.array([float(r["relative_l2_error"]) * 100 for r in points])
            bars = np.array([[float(r["loop_median_ms"]) - float(r["loop_q25_ms"]) for r in points],
                             [float(r["loop_q75_ms"]) - float(r["loop_median_ms"]) for r in points]])
            ax.errorbar(x, y, xerr=bars, fmt=marker, color=color, label=label, capsize=3)
            for r, a, b in zip(points, x, y):
                ax.annotate(f"{float(r['dt_ms']):g} ms", (a, b), xytext=(5, 5),
                            textcoords="offset points", fontsize=8)
        ax.set(xlabel="Boucle médiane (ms)", ylabel="Erreur finale relative L² (%)",
               title=title, xscale="log", yscale="log")
        ax.grid(alpha=.2)
        ax.legend()
    for budget, color in [(0.5, "#835caa"), (1., "#b96d2a"), (2., "#078578"), (5., "#286b9c")]:
        points = sorted([r for r in frontier if r["budget_percent"] == budget
                         and r["scenario"] != "original_three_layers"], key=lambda r: r["fine_fraction"])
        axes[2].plot([r["fine_fraction"] * 100 for r in points],
                     [r["speedup_global_over_local"] for r in points], "o-", color=color, label=f"Erreur ≤ {budget:g} %")
    axes[2].axhline(1., color="black", lw=1, ls="--")
    if refined:
        points = sorted([r for r in refined if r["scenario"] != "original_three_layers"],
                        key=lambda r: r["fine_fraction"])
        axes[2].plot([r["fine_fraction"] * 100 for r in points],
                     [r["speedup_global_over_local"] for r in points], "x--", color="black", label="2 %, pas raffiné")
    axes[2].set(xlabel="Éléments au petit pas (%)", ylabel="Temps global / temps local",
                title="Plus rapide parmi les pas testés\nAu-dessus de 1 : avantage au local")
    axes[2].grid(alpha=.2)
    axes[2].legend(fontsize=9)
    fig.supxlabel(f"Kokkos {protocol['backend']} · {protocol['repetitions']} répétitions · 1 cœur fixé · Durée physique 1,8 s\n"
                  "Référence globale : 0,375 ms sur le même maillage · Barres : intervalle interquartile · Initialisation exclue",
                  fontsize=9)
    fig.savefig(root / "06_fraction_precision_performance.png", dpi=160)
    print(f"Validated {len(samples)} samples, {len(summary)} configurations and {len(frontier)} budget comparisons.")
    print("scenario | fine % | best local at 2% | best global at 2% | global/local")
    for r in frontier:
        if r["budget_percent"] == 2.:
            print(f"{r['scenario']} | {r['fine_fraction'] * 100:.1f} | {r['local_case']} {r['local_ms']:.2f} ms"
                  f" | {r['global_case']} {r['global_ms']:.2f} ms | {r['speedup_global_over_local']:.3f}")


if __name__ == "__main__":
    main()
