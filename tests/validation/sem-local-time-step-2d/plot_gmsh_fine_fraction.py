"""Figure and table of the fine-fraction study (sem2d_gmsh_fine_fraction): gain of the
two-level local time step over the global scheme at its stability limit, against the
fraction of fine elements, on CPU and GPU, with the bound 2 / (1 + support fraction)."""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RESULTS = Path(__file__).resolve().parents[3] / "documentation" / "acoustic" / "local-time-step"


def read(path):
    with open(path, newline="") as f:
        return [{k: (v if k == "case" else float(v) if v != "" else float("nan")) for k, v in row.items()}
                for row in csv.DictReader(f)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS / "fine-fraction-gmsh-results")
    args = parser.parse_args()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    lines = ["| Backend | Cas | Éléments | Fins | Support | Global (ms) | LTS (ms) | Gain | Gain, même noyau | Borne "
             "| Erreur globale | Erreur LTS |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    markers = {"40": "o", "20": "s", "10": "^"}
    for backend, color in (("cpu", "C0"), ("gpu", "C3")):
        path = args.results / f"{backend}.csv"
        if not path.exists():
            continue
        rows = read(path)
        for h, marker in markers.items():
            # side0 = no patch: every element is fine (degenerate), plotted in the table only
            sel = sorted((r for r in rows if r["case"].startswith(f"h{h}_") and not r["case"].endswith("side0")
                          and r["stable"] == 1), key=lambda r: r["fine_fraction"])
            if not sel:
                continue
            x = [100 * r["fine_fraction"] for r in sel]
            gain = [r["gain"] for r in sel]
            # interquartile range of the gain: (global q25 / local q75, global q75 / local q25)
            low = [g - r["global_ms_q25"] / r["local_ms_q75"] for g, r in zip(gain, sel)]
            high = [r["global_ms_q75"] / r["local_ms_q25"] - g for g, r in zip(gain, sel)]
            ls = "-" if backend == "cpu" else "none"   # GPU timings too dispersed to draw a trend
            ax.errorbar(x, gain, yerr=[low, high], marker=marker, color=color, ls=ls, capsize=3,
                        label=f"{backend.upper()}, h = {h} m")
            if backend == "cpu":
                ax.plot(x, [r["gain_same_kernel"] for r in sel], marker=marker, color=color, ls=":", mfc="none",
                        label=f"CPU, h = {h} m, noyau de correction de sem2d")
        for r in rows:
            bound = 2 / (1 + r["support_fraction"])
            if r["stable"] != 1:
                lines.append(f"| {backend.upper()} | {r['case']} | {int(r['elements'])} | {100 * r['fine_fraction']:.1f} % | "
                             f"{100 * r['support_fraction']:.1f} % | LTS instable | | | | {bound:.2f} | | |")
                continue
            lines.append(f"| {backend.upper()} | {r['case']} | {int(r['elements'])} | {100 * r['fine_fraction']:.1f} % | "
                         f"{100 * r['support_fraction']:.1f} % | {r['global_ms']:.0f} | {r['local_ms']:.0f} | "
                         f"{r['gain']:.2f} | {r['gain_same_kernel']:.2f} | {bound:.2f} | {r['error_global']:.1e} | "
                         f"{r['error_local']:.1e} |")
    import numpy as np
    s = np.linspace(0, 0.5, 100)
    ax.plot(100 * s, 2 / (1 + s), "k--", lw=1, label="borne 2 / (1 + support), support ≈ fins")
    ax.axhline(1, color="0.5", lw=0.8)
    ax.set_xlabel("Éléments au petit pas (%)")
    ax.set_ylabel("Gain (temps global / temps LTS)")
    ax.set_title("Pas de temps local à deux niveaux, maillages gmsh raffinés")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(args.results / "gain_vs_fine_fraction.png", dpi=150)
    (args.results / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
