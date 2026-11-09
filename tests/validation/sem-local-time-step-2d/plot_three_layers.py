"""Render measured three-layer SEM2D LTS validation data; no simulation in Python."""

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/funtides-lts-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import MaxNLocator
import numpy as np

# Results and reports of the validations
RESULTS = Path(__file__).resolve().parents[3] / "documentation" / "acoustic" / "local-time-step"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS / "results")
    args = parser.parse_args()
    root = args.results
    meta = json.loads((root / "metadata.json").read_text())
    elements = np.genfromtxt(root / "elements.csv", delimiter=",", names=True)
    fields = np.genfromtxt(root / "fields.csv", delimiter=",", names=True)
    history = np.genfromtxt(root / "history.csv", delimiter=",", names=True)
    assert np.isfinite(fields["local"]).all() and np.isfinite(fields["reference"]).all()
    assert np.isfinite(history["energy"]).all()
    nx, nz = meta["nx"], meta["nz"]
    frames = [fields[fields["frame"] == frame] for frame in np.unique(fields["frame"])]
    assert all(len(frame) == nx * nz for frame in frames)
    x = frames[0]["x"].reshape(nz, nx)[0] / 1000
    z = frames[0]["z"].reshape(nz, nx)[:, 0] / 1000

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10.5,
        "axes.titlesize": 12, "axes.labelsize": 10.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })

    def save(fig, name):
        fig.savefig(root / name, dpi=160, bbox_inches="tight")
        plt.close(fig)
        print(root / name)

    def interfaces(ax, labels=False):
        for height in (1, 2):
            ax.axhline(height, color="#172c42", lw=0.8, ls="--", alpha=0.55)
        ax.set(xlim=(0, 6), ylim=(0, 3), xlabel="x (km)", ylabel="z (km, vers le haut)")
        ax.set_aspect("equal")

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), layout="constrained")
    fig.suptitle("Trois couches : vitesse et pas de temps local", fontsize=16, fontweight="bold", y=1.08)
    vp = elements["vp"].reshape(30, 60)
    cm = ListedColormap(["#264c7b", "#279eaa", "#e6c650"])
    im = axes[0].imshow(vp, extent=(0, 6, 0, 3), origin="lower", interpolation="nearest",
                        cmap=cm, norm=BoundaryNorm([1125, 1875, 2625, 3375], 3), aspect="equal")
    axes[0].set_title("Même maille de 100 m dans les trois couches")
    interfaces(axes[0])
    for height, speed in [(2.5, 1500), (1.5, 2250), (0.5, 3000)]:
        axes[0].text(4.5, height, f"{speed:,} m/s".replace(",", " "), ha="center", va="center",
                     bbox=dict(facecolor="white", alpha=0.88, edgecolor="none", pad=5))
    axes[0].plot(3, 2.5, marker="*", color="#f73551", markersize=15, markeredgecolor="white")
    fig.colorbar(im, ax=axes[0], orientation="horizontal", shrink=0.85, pad=0.04,
                 ticks=[1500, 2250, 3000], label="Vitesse de propagation Vp (m/s)")
    steps = elements["selected_dt"].reshape(30, 60) * 1000
    dt_cm = ListedColormap(["#80cbe5", "#f1b95b"])
    im = axes[1].imshow(steps, extent=(0, 6, 0, 3), origin="lower", interpolation="nearest",
                        cmap=dt_cm, norm=BoundaryNorm([3, 9, 15], 2), aspect="equal")
    axes[1].set_title("Deux niveaux temporels, rapport 2")
    interfaces(axes[1])
    for height, estimate, selected in [(2.5, 13.333, 12), (1.5, 8.889, 6), (0.5, 6.667, 6)]:
        axes[1].text(3, height, f"Estimé : {estimate:.2f} ms   →   Retenu : {selected} ms",
                     ha="center", va="center", bbox=dict(facecolor="white", alpha=0.9, edgecolor="none", pad=5))
    fig.colorbar(im, ax=axes[1], orientation="horizontal", shrink=0.85, pad=0.04, ticks=[6, 12],
                 label="Pas attribué aux éléments (ms)")
    fig.supxlabel("SEM ordre 2 · CFL de partition 0,4 · Les nœuds partagés suivent la zone à petit pas.", fontsize=10)
    save(fig, "01_modele_et_pas_locaux.png")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8.1), layout="constrained")
    fig.suptitle("Propagation calculée avec le pas de temps local", fontsize=16, fontweight="bold")
    peak = max(np.abs(frame["local"]).max() for frame in frames[1:])
    limit = np.ceil(peak * 100) / 100
    for ax, frame in zip(axes.flat, frames[1:]):
        pressure = frame["local"].reshape(nz, nx)
        im = ax.pcolormesh(x, z, pressure, shading="gouraud", cmap="RdBu_r", vmin=-limit, vmax=limit,
                           rasterized=True)
        interfaces(ax)
        ax.set_title(f"t = {frame['time'][0]:.2f} s")
    fig.colorbar(im, ax=list(axes.flat), shrink=0.85, pad=0.02, label="Pression (unités arbitraires)")
    fig.supxlabel("Impulsion initiale gaussienne à (3 km ; 2,5 km), σ = 160 m · Bords rigides · Sans source ni amortissement",
                  fontsize=9.5)
    save(fig, "02_propagation_trois_couches.png")

    final = frames[-1]
    a = final["local"].reshape(nz, nx)
    b = final["reference"].reshape(nz, nx)
    difference = a - b
    common = max(np.abs(a).max(), np.abs(b).max())
    error_limit = np.abs(difference).max()
    fig, axes = plt.subplots(1, 3, figsize=(15.3, 4.6), layout="constrained")
    fig.suptitle("Validation à t = 1,80 s : pas local contre petit pas global", fontsize=15, fontweight="bold")
    for ax, values, title in zip(axes, [b, a, difference], ["Référence globale · 1,5 ms", "Pas local · 12 / 6 ms", "Différence : local − référence"]):
        lim = error_limit if ax is axes[2] else common
        im = ax.pcolormesh(x, z, values, shading="gouraud", cmap="RdBu_r", vmin=-lim, vmax=lim,
                           rasterized=True)
        interfaces(ax)
        ax.set_title(title)
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.05, shrink=0.9,
                          label="Échelle propre à l’écart (u.a.)" if ax is axes[2] else "Pression (u.a.)")
        cb.locator = MaxNLocator(nbins=5)
        cb.update_ticks()
    fig.supxlabel(f"Erreur relative L² pondérée par la masse : {meta['relative_l2_error'] * 100:.2f} % · Même maillage pour les deux calculs",
                  fontsize=11)
    save(fig, "03_comparaison_reference.png")

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.5), layout="constrained")
    fig.suptitle("Contrôles quantitatifs du cas à trois couches", fontsize=15, fontweight="bold")
    axes[0].plot(history["time"], 100 * history["relative_l2_error"], color="#266899", lw=2)
    axes[0].set(xlabel="Temps (s)", ylabel="Erreur relative L² pondérée par la masse (%)",
                title="Écart au calcul global à 1,5 ms", xlim=(0, history["time"][-1]))
    axes[0].annotate(f"{meta['relative_l2_error'] * 100:.2f} % à 1,80 s",
                     xy=(history["time"][-1], 100 * history["relative_l2_error"][-1]),
                     xytext=(-12, -35), textcoords="offset points", ha="right",
                     arrowprops=dict(arrowstyle="->", color="#266899"))
    drift = history["relative_energy_drift"] * 1e6
    axes[1].plot(history["time"], drift, color="#078578", lw=1.6)
    axes[1].axhline(0, color="#65707c", lw=0.8)
    axes[1].set(xlabel="Temps (s)", ylabel="(E / E₀ − 1) × 10⁶",
                title=f"Énergie discrète · dérive max. {meta['max_relative_energy_drift']:.2e}",
                xlim=(0, history["time"][-1]))
    for ax in axes:
        ax.grid(alpha=0.22)
        ax.set_axisbelow(True)
    fig.supxlabel("CPU Kokkos Serial · 1 800 éléments · 7 381 nœuds · 150 grands pas locaux / 1 200 pas de référence",
                  fontsize=10)
    save(fig, "04_erreur_et_energie.png")


if __name__ == "__main__":
    main()
