"""Plot repeated compute-time measurements and their accuracy trade-off."""

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/funtides-lts-matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path(__file__).parent / "results")
    args = parser.parse_args()
    rows = np.genfromtxt(args.results / "timing_summary.csv", delimiter=",", names=True, dtype=None, encoding="utf8")
    protocol = json.loads((args.results / "timing_protocol.json").read_text())
    names = ["Local 12 / 6 ms", "Global 9 ms", "Global 6 ms", "Global 3 ms", "Global 1,5 ms"]
    colors = ["#078578", "#286b9c", "#7695ad", "#94a9ba", "#b3bfc9"]
    values = rows["loop_median_ms"]
    errors = rows["relative_l2_error"] * 100
    whiskers = np.vstack((values - rows["loop_q25_ms"], rows["loop_q75_ms"] - values))
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    fig.suptitle("Cas à trois couches : temps de calcul et précision", fontsize=16, fontweight="bold")
    yy = np.arange(len(rows))
    axes[0].barh(yy, values, color=colors, xerr=whiskers, error_kw={"capsize": 3, "lw": 1})
    axes[0].set_yticks(yy, names)
    axes[0].invert_yaxis()
    axes[0].set(xlabel="Temps médian de la boucle (ms)", title="Même durée physique : 1,8 s",
                xlim=(0, values.max() * 1.2))
    for index, value in enumerate(values):
        axes[0].text(value + values.max() * 0.025, index, f"{value:.1f} ms", va="center", fontsize=10)
    offsets = [(9, 7), (8, -17), (8, 7), (8, 7), (-8, 8)]
    for i in range(len(rows)):
        axes[1].scatter(values[i], errors[i], s=85, color=colors[i], edgecolor="white", zorder=3)
        axes[1].annotate(names[i], (values[i], errors[i]), xytext=offsets[i], textcoords="offset points",
                         ha="right" if i == 4 else "left", fontsize=10)
    axes[1].set(xlabel="Temps médian de la boucle (ms)", ylabel="Erreur relative L² (%)",
                title="Erreur face au calcul global à 0,75 ms", xscale="log", yscale="log")
    axes[1].set_xlim(values.min() * .75, values.max() * 1.35)
    axes[1].set_ylim(errors.min() * .6, errors.max() * 1.8)
    for ax in axes:
        ax.grid(alpha=.2, axis="x" if ax is axes[0] else "both")
        ax.set_axisbelow(True)
    fig.supxlabel(f"CPU Kokkos {protocol['backend']} · 1 cœur fixé · {protocol['repetitions']} répétitions après échauffement · Barres : intervalle interquartile\n"
                  "Chronométrage sans initialisation, calcul d’énergie, sorties ni copies vers le host.", fontsize=9.5)
    path = args.results / "05_temps_calcul.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    print(path)


if __name__ == "__main__":
    main()
