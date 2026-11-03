#!/usr/bin/env python3
"""Resume productos IRAF por chip y crea mapas y diagnosticos de fotometria.

Requiere numpy, matplotlib y astropy. Ejecutar desde cualquier carpeta:
    python3 analisis_fotometria.py [--max-error-percent 10]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval
from scipy.spatial import cKDTree


def read_quality(path: Path) -> list[dict[str, float]]:
    """Lee la tabla x y msky stdev nsky flux mag merr (con o sin cabecera)."""
    rows = []
    with path.open() as stream:
        for line in stream:
            bits = line.split()
            if not bits or bits[0].startswith("#") or bits[0].lower() == "x":
                continue
            if len(bits) < 8:
                continue
            try:
                x, y, sky, stdev, nsky, flux, mag, merr = map(float, bits[:8])
            except ValueError:
                continue
            relerr = 100.0 * (math.log(10.0) / 2.5) * merr if math.isfinite(merr) else math.inf
            measured = math.isfinite(flux) and flux > 0 and math.isfinite(mag)
            rows.append({"x": x, "y": y, "msky": sky, "stdev": stdev,
                         "nsky": nsky, "flux": flux, "mag": mag,
                         "merr": merr, "error_percent": relerr,
                         "measured": measured})
    return rows


def image_data(path: Path) -> np.ndarray:
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            if hdu.data is not None and np.asarray(hdu.data).ndim == 2:
                return np.asarray(hdu.data, dtype=np.float32)
    raise ValueError(f"No se encontro imagen 2D en {path}")


def nearest_distances(rows: list[dict[str, float]], indexes: list[int]) -> dict[int, float]:
    if len(indexes) < 2:
        return {idx: float("nan") for idx in indexes}
    xy = np.array([[rows[i]["x"], rows[i]["y"]] for i in indexes])
    distances, _ = cKDTree(xy).query(xy, k=2)
    return {idx: float(distances[j, 1]) for j, idx in enumerate(indexes)}


def make_chip_outputs(chipdir: Path, max_error: float) -> dict:
    tag = chipdir.name
    prefix = f"V_{tag}"
    catalog = chipdir / f"{prefix}_quality.dat"
    fits_path = chipdir / f"SexA_V_{tag}.fits"
    if not catalog.exists() or not fits_path.exists():
        return {"chip": tag, "estado": "omitido: falta quality.dat o imagen FITS"}

    rows = read_quality(catalog)
    valid_idx = [i for i, row in enumerate(rows) if row["measured"]]
    approved_idx = [i for i in valid_idx if rows[i]["error_percent"] <= max_error]
    distances = nearest_distances(rows, approved_idx)
    for idx, row in enumerate(rows):
        row["vecino_mas_cercano_px"] = distances.get(idx, float("nan"))
        row["aprobada"] = idx in approved_idx

    csv_path = chipdir / f"{prefix}_analisis_porcentual.csv"
    columns = ["x", "y", "flux", "mag", "merr", "error_percent", "vecino_mas_cercano_px", "aprobada"]
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in columns})

    img = image_data(fits_path)
    norm = ImageNormalize(img, interval=ZScaleInterval(), stretch=AsinhStretch())
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    ax.imshow(img, origin="lower", cmap="gray", extent=(1, img.shape[1], 1, img.shape[0]),
              norm=norm, interpolation="nearest")
    if approved_idx:
        xy = np.array([[rows[i]["x"], rows[i]["y"]] for i in approved_idx])
        ax.scatter(xy[:, 0], xy[:, 1], s=30, facecolors="none", edgecolors="#00ff55",
                   linewidths=0.9, label=f"Aprobadas ({len(approved_idx)})")
        ax.legend(loc="upper right", facecolor="black", labelcolor="white", framealpha=0.7)
    ax.set(title=f"{tag}: fotometría aprobada (error relativo ≤ {max_error:g}%)",
           xlabel="X [pixeles]", ylabel="Y [pixeles]")
    fig.savefig(chipdir / f"{prefix}_estrellas_aprobadas.png", dpi=160)
    plt.close(fig)

    nn = np.array([distances[i] for i in approved_idx if math.isfinite(distances[i])])
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    if nn.size:
        ax.hist(nn, bins="auto", color="steelblue", edgecolor="white")
        ax.axvline(float(np.median(nn)), color="darkorange", linestyle="--",
                   label=f"Mediana: {np.median(nn):.2f} pix")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "No hay suficientes estrellas aprobadas para distancias",
                ha="center", va="center", transform=ax.transAxes)
    ax.set(title=f"{tag}: distancia al vecino aprobado más cercano",
           xlabel="Distancia al vecino más cercano [pixeles]", ylabel="Número de estrellas")
    fig.savefig(chipdir / f"{prefix}_histograma_vecinos.png", dpi=160)
    plt.close(fig)

    summary = {
        "chip": tag,
        "catalogo": catalog.name,
        "imagen": fits_path.name,
        "error_maximo_porcentual": max_error,
        "estrellas_con_fotometria_valida": len(valid_idx),
        "estrellas_aprobadas": len(approved_idx),
        "estrellas_rechazadas_por_error": len(valid_idx) - len(approved_idx),
        "error_relativo_mediano_aprobadas_pct": float(np.median([rows[i]["error_percent"] for i in approved_idx])) if approved_idx else None,
        "vecino_mas_cercano_mediano_pix": float(np.median(nn)) if nn.size else None,
        "tabla": csv_path.name,
        "imagen_estrellas": f"{prefix}_estrellas_aprobadas.png",
        "histograma_vecinos": f"{prefix}_histograma_vecinos.png",
    }
    (chipdir / f"{prefix}_resumen_fotometria.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-error-percent", type=float, default=10.0,
                        help="máximo error relativo de flujo para aprobar una estrella (default: 10%%)")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent,
                        help="directorio del proyecto que contiene chips/")
    args = parser.parse_args()
    if args.max_error_percent < 0:
        parser.error("--max-error-percent debe ser no negativo")
    summaries = [make_chip_outputs(p, args.max_error_percent)
                 for p in sorted((args.root / "chips").glob("chip*")) if p.is_dir()]
    completed = [s for s in summaries if "estrellas_con_fotometria_valida" in s]
    total_measured = sum(s["estrellas_con_fotometria_valida"] for s in completed)
    total_approved = sum(s["estrellas_aprobadas"] for s in completed)
    report = {"error_maximo_porcentual": args.max_error_percent,
              "total_estrellas_con_fotometria_valida": total_measured,
              "total_estrellas_aprobadas": total_approved,
              "chips": summaries}
    out = args.root / "chips" / "resumen_fotometria_porcentual.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"Fotometría válida: {total_measured} estrellas")
    print(f"Aprobadas con error ≤ {args.max_error_percent:g}%: {total_approved} estrellas")
    print(f"Resumen: {out}")


if __name__ == "__main__":
    main()
