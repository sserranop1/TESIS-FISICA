#!/usr/bin/env python3
"""Regenerate the per-filter chip 5 comparison figures for the thesis."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize


ROOT = Path(__file__).resolve().parents[1]
CHIP_DIR = ROOT / "dataSextansA_multiplicadores_FWHM" / "chips" / "chip5"
IMAGE_DIR = ROOT / "dataSextansA" / "chips" / "chip5"
OUTPUT_DIR = Path(__file__).resolve().parent / "img"


def read_xy(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2:
            try:
                rows.append((float(fields[0]), float(fields[1])))
            except ValueError:
                continue
    return np.asarray(rows, dtype=float).reshape((-1, 2))


def read_valid_photometry(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 8:
            continue
        try:
            x, y = float(fields[0]), float(fields[1])
            mag, merr = float(fields[6]), float(fields[7])
        except ValueError:
            continue
        if np.isfinite(mag) and np.isfinite(merr):
            rows.append((x, y))
    return np.asarray(rows, dtype=float).reshape((-1, 2))


def read_traditional_xy(chip_dir: Path, filt: str) -> np.ndarray:
    groups = chip_dir / f"{filt}_chip5_grupos_adaptativos.csv"
    with groups.open(newline="") as stream:
        first_group = next(csv.DictReader(stream))
    # The first group is the unchanged, initially accepted fixed annulus.
    return read_xy(chip_dir / Path(first_group["coords"]).name)


def read_recovered_xy(chip_dir: Path, filt: str) -> np.ndarray:
    path = chip_dir / f"{filt}_chip5_contaminacion.csv"
    rows = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            rows.append((float(row["x"]), float(row["y"])))
    return np.asarray(rows, dtype=float).reshape((-1, 2))


def load_image(filt: str) -> np.ndarray:
    path = IMAGE_DIR / f"SexA_{filt}_chip5.fits"
    with fits.open(path, memmap=True) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float32).copy()
    data[~np.isfinite(data) | (data >= 49999)] = np.nan
    return data


def image_norm(data: np.ndarray) -> ImageNormalize:
    valid = data[np.isfinite(data)]
    low, high = np.percentile(valid, [0.5, 99.6])
    return ImageNormalize(vmin=low, vmax=high, stretch=AsinhStretch(a=0.03), clip=True)


def save_overview(filt: str, data: np.ndarray, phot: np.ndarray,
                  traditional: np.ndarray, recovered: np.ndarray) -> None:
    height, width = data.shape
    cmap = plt.get_cmap("gray").copy()
    cmap.set_bad("#080d14")
    fig, ax = plt.subplots(figsize=(16, 8.15), constrained_layout=True)
    ax.imshow(data, origin="lower", cmap=cmap, norm=image_norm(data),
              extent=(0.5, width + 0.5, 0.5, height + 0.5), interpolation="nearest")

    if len(phot):
        ax.scatter(phot[:, 0], phot[:, 1], s=7, marker=".", color="#32c9ef",
                   alpha=0.55, linewidths=0, rasterized=True,
                   label=f"Fotometría válida ({len(phot):,})")
    if len(traditional):
        ax.scatter(traditional[:, 0], traditional[:, 1], s=22, marker="o",
                   facecolors="none", edgecolors="#62f06a", linewidths=0.65,
                   rasterized=True, label=f"Anillo fijo aceptado ({len(traditional):,})")
    if len(recovered):
        ax.scatter(recovered[:, 0], recovered[:, 1], s=78, marker="*",
                   facecolors="#ff9d22", edgecolors="#161616", linewidths=0.7,
                   zorder=5, label=f"Recuperadas con anillo adaptativo ({len(recovered):,})")

    ax.set(xlim=(0, width), ylim=(0, height), xlabel="x [píxeles]", ylabel="y [píxeles]",
           title=f"Sextans A · chip 5 · filtro {filt}")
    ax.legend(loc="upper right", framealpha=0.88, fontsize=12, markerscale=1.2)
    output = OUTPUT_DIR / f"sextansA_chip5_clasificacion_{filt}.png"
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def read_moved_example(chip_dir: Path, filt: str) -> dict[str, str]:
    with (chip_dir / f"{filt}_chip5_movidas.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    # Use the V example already discussed in the text and the previous figure.
    return rows[1] if filt == "V" else rows[0]


def save_cutout(filt: str, data: np.ndarray, row: dict[str, str]) -> None:
    x, y = float(row["x"]), float(row["y"])
    ap_radius = {"B": 9.26, "V": 8.96, "I": 8.22}[filt]
    fixed_inner = {"B": 23.15, "V": 22.40, "I": 20.55}[filt]
    fixed_outer = {"B": 38.15, "V": 37.40, "I": 35.55}[filt]
    adaptive_inner, adaptive_outer = float(row["r_in"]), float(row["r_out"])
    half = int(np.ceil(adaptive_outer + 23))

    height, width = data.shape
    x0, x1 = max(0, int(x - half)), min(width, int(x + half))
    y0, y1 = max(0, int(y - half)), min(height, int(y + half))
    cut = data[y0:y1, x0:x1]
    cmap = plt.get_cmap("gray").copy()
    cmap.set_bad("#080d14")

    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    ax.imshow(cut, origin="lower", cmap=cmap, norm=image_norm(cut),
              extent=(x0 + 0.5, x1 + 0.5, y0 + 0.5, y1 + 0.5), interpolation="nearest")
    ax.add_patch(plt.Circle((x, y), ap_radius, fill=False, color="#29d9f2",
                            linewidth=2.0, linestyle=":"))
    for radius in (fixed_inner, fixed_outer):
        ax.add_patch(plt.Circle((x, y), radius, fill=False, color="#ffad28",
                                linewidth=1.8, linestyle="--"))
    for radius in (adaptive_inner, adaptive_outer):
        ax.add_patch(plt.Circle((x, y), radius, fill=False, color="#58ef76",
                                linewidth=1.8, linestyle="-"))
    ax.plot(x, y, marker="+", markersize=12, markeredgewidth=2, color="white")

    from matplotlib.lines import Line2D

    legend = [
        Line2D([0], [0], color="#29d9f2", linestyle=":", lw=2, label="Apertura"),
        Line2D([0], [0], color="#ffad28", linestyle="--", lw=2, label="Anillo fijo"),
        Line2D([0], [0], color="#58ef76", linestyle="-", lw=2, label="Anillo adaptativo"),
    ]
    ax.legend(handles=legend, loc="upper right", framealpha=0.9)
    ax.set(xlim=(x - half, x + half), ylim=(y - half, y + half),
           xlabel="x [píxeles]", ylabel="y [píxeles]",
           title=(f"Filtro {filt} · fuente recuperada · "
                  f"$\\Delta MAG={float(row['delta_mag']):+.3f}$ mag"))
    output = OUTPUT_DIR / f"sextansA_chip5_anillos_{filt}.png"
    fig.savefig(output, dpi=260, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for filt in "BVI":
        phot = read_valid_photometry(CHIP_DIR / f"{filt}_chip5_fixed.dat")
        traditional = read_traditional_xy(CHIP_DIR, filt)
        recovered = read_recovered_xy(CHIP_DIR, filt)
        data = load_image(filt)
        save_overview(filt, data, phot, traditional, recovered)
        save_cutout(filt, data, read_moved_example(CHIP_DIR, filt))
        print(f"{filt}: fotometría={len(phot)}, tradicional={len(traditional)}, recuperadas={len(recovered)}")


if __name__ == "__main__":
    main()
