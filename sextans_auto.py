#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sextans A | Automatización completa con IRAF/DAOFIND/PHOT.
Derivado de las funciones del cuaderno de tesis suministrado el 26-10-2026.
Selección automática de 15 estrellas y 15 cajas de cielo distribuidas.
Usar --help. Requiere astropy, scipy, pandas, numpy, matplotlib, IRAF.
"""
import argparse
import json
import logging
import math
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.spatial import cKDTree
from scipy.optimize import least_squares



# Ruta principal del proyecto en el Mac
BASE = (
    Path.home()
    / "Desktop"
    / "Universidad"
    / "tesis"
    / "dataSextansA"
)

THRESHOLDS_DEFAULT = np.array([
    2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0
])

MERR_MAX_DEFAULT = 0.20
ANCHO_SKY_DEFAULT = 15.0
PASO_ANILLO_DEFAULT = 5.0
N_ANILLOS_DEFAULT = 5
VALOR_INVALIDO_DEFAULT = 50000.0


# === FUNCIONES DEL NOTEBOOK, CELDA 110 ===
def resumir_fwhm(valores, mostrar=True):
    """
    Calcula los estadísticos del FWHM a partir de las estrellas
    que ya fueron seleccionadas manualmente como adecuadas.
    """
    valores = np.asarray(valores, dtype=float)
    valores = valores[np.isfinite(valores)]

    if len(valores) == 0:
        raise ValueError("No hay valores válidos de FWHM.")

    resultado = {
        "n": len(valores),
        "mediana": float(np.median(valores)),
        "media": float(np.mean(valores)),
        "std": float(np.std(valores, ddof=1)) if len(valores) > 1 else np.nan,
    }

    resultado["r_ap"] = 3.0 * resultado["mediana"]
    resultado["r_sky_in"] = 5.0 * resultado["mediana"]
    resultado["r_sky_out"] = (
        resultado["r_sky_in"] + ANCHO_SKY_DEFAULT
    )

    if mostrar:
        print(f"N estrellas        = {resultado['n']}")
        print(f"FWHM mediano       = {resultado['mediana']:.4f} pix")
        print(f"FWHM medio         = {resultado['media']:.4f} pix")
        print(f"Desv. estándar     = {resultado['std']:.4f} pix")
        print(f"Radio apertura     = {resultado['r_ap']:.4f} pix")
        print(f"Radio sky interno  = {resultado['r_sky_in']:.4f} pix")
        print(f"Radio sky externo  = {resultado['r_sky_out']:.4f} pix")

    return resultado


def resumir_cielo(sigmas, medianas=None, mostrar=True):
    """
    Aplica Tukey a la distribución de sigma de las cajas de cielo.
    Devuelve sigma_sky como la mediana de las cajas aceptadas.
    """
    sigmas = np.asarray(sigmas, dtype=float)
    sigmas = sigmas[np.isfinite(sigmas)]

    if len(sigmas) == 0:
        raise ValueError("No hay mediciones válidas de sigma del cielo.")

    q1 = np.percentile(sigmas, 25)
    q3 = np.percentile(sigmas, 75)
    iqr = q3 - q1
    limite = q3 + 1.5 * iqr

    mask = sigmas <= limite
    aceptadas = sigmas[mask]
    rechazadas = sigmas[~mask]

    resultado = {
        "q1": float(q1),
        "q3": float(q3),
        "iqr": float(iqr),
        "limite": float(limite),
        "n_total": len(sigmas),
        "n_aceptadas": len(aceptadas),
        "n_rechazadas": len(rechazadas),
        "sigma_sky": float(np.median(aceptadas)),
        "aceptadas": aceptadas,
        "rechazadas": rechazadas,
    }

    if medianas is not None:
        med = np.asarray(medianas, dtype=float)
        med = med[np.isfinite(med)]
        if len(med):
            resultado["sky_mediana"] = float(np.median(med))

    if mostrar:
        print(f"Q1                 = {q1:.4f} ADU")
        print(f"Q3                 = {q3:.4f} ADU")
        print(f"IQR                = {iqr:.4f} ADU")
        print(f"Límite Tukey       = {limite:.4f} ADU")
        print(f"Cajas aceptadas    = {len(aceptadas)}")
        print(f"Cajas rechazadas   = {len(rechazadas)}")
        print(f"Sigma sky adoptado = {resultado['sigma_sky']:.4f} ADU")
        if "sky_mediana" in resultado:
            print(f"Mediana del cielo  = {resultado['sky_mediana']:.4f} ADU")

    return resultado


# === FUNCIONES DEL NOTEBOOK, CELDA 112 ===
def rutas_chip(chip, filtro="V", base=BASE):
    chip = int(chip)
    carpeta = base / "chips" / f"chip{chip}"

    return {
        "carpeta": carpeta,
        "imagen": carpeta / f"SexA_{filtro}_chip{chip}.fits",
        "prefijo": f"{filtro}_chip{chip}",
    }


def generar_script_daofind(
    chip,
    filtro,
    fwhm,
    sigma_sky,
    thresholds=THRESHOLDS_DEFAULT,
    base=BASE,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    imagen = rutas["imagen"]
    prefijo = rutas["prefijo"]

    hdr = fits.getheader(imagen)

    gain = hdr.get("GAIN", 2.75)
    rdnoise = hdr.get("RDNOISE", 6.4125)
    saturate = hdr.get("SATURATE", 39000)

    lineas = [
        f'datapars.fwhmpsf = {fwhm:.6f}',
        f'datapars.sigma = {sigma_sky:.6f}',
        'datapars.datamin = INDEF',
        f'datapars.datamax = {saturate}',
        f'datapars.readnoise = {rdnoise}',
        f'datapars.epadu = {gain}',
        '',
    ]

    for th in thresholds:
        tag = str(th).replace(".", "_")
        salida = carpeta / f"{prefijo}_thr{tag}.coo"

        lineas += [
            f'findpars.threshold = {th}',
            (
                f'daofind {imagen} '
                f'output={salida} verify-'
            ),
            '',
        ]

    script = carpeta / f"{prefijo}_daofind.cl"
    script.write_text("\n".join(lineas), encoding="utf-8")

    print("Script creado:")
    print(script)
    print()
    print("Ejecutar dentro de IRAF:")
    print(f"cl < {script}")

    return script


# === FUNCIONES DEL NOTEBOOK, CELDA 114 ===
def analizar_threshold(thresholds, detecciones, graficar=True):
    thresholds = np.asarray(thresholds, dtype=float)
    detecciones = np.asarray(detecciones, dtype=float)

    if len(thresholds) != len(detecciones):
        raise ValueError("thresholds y detecciones deben tener igual longitud.")

    if len(thresholds) < 3 or not np.isfinite(thresholds).all():
        raise ValueError("Se requieren al menos tres thresholds finitos.")

    if len(np.unique(thresholds)) < 3:
        raise ValueError("Se requieren al menos tres thresholds distintos.")

    if not np.isfinite(detecciones).all():
        raise ValueError("Todos los conteos deben ser finitos.")

    if np.any(detecciones <= 0):
        raise ValueError("Todos los conteos deben ser positivos.")

    lnN = np.log(detecciones)

    a, b, c = np.polyfit(thresholds, lnN, 2)
    if not np.isfinite([a, b, c]).all() or np.isclose(a, 0.0):
        raise ValueError(
            "El ajuste cuadrático no define un threshold óptimo finito."
        )

    ajuste = np.polyval([a, b, c], thresholds)

    ss_res = np.sum((lnN - ajuste) ** 2)
    ss_tot = np.sum((lnN - np.mean(lnN)) ** 2)
    if np.isclose(ss_tot, 0.0):
        raise ValueError("Los conteos no varían entre thresholds.")
    r2 = 1 - ss_res / ss_tot

    threshold_opt = (-1 - b) / (2 * a)
    if not np.isfinite(threshold_opt):
        raise ValueError("El threshold óptimo calculado no es finito.")

    tabla = pd.DataFrame({
        "threshold": thresholds,
        "detecciones": detecciones.astype(int),
        "lnN": lnN,
    })

    resultado = {
        "a": float(a),
        "b": float(b),
        "c": float(c),
        "r2": float(r2),
        "threshold_optimo": float(threshold_opt),
        "tabla": tabla,
    }

    print(
        f"ln(N) = {a:.8f} x² "
        f"{b:+.8f} x {c:+.8f}"
    )
    print(f"R² = {r2:.8f}")
    print(f"Threshold calculado = {threshold_opt:.4f}")

    if graficar:
        plt.figure(figsize=(7, 4.5))
        plt.plot(thresholds, detecciones, marker="o")
        plt.axvline(
            threshold_opt,
            linestyle="--",
            label=f"Threshold = {threshold_opt:.2f}"
        )
        plt.xlabel("Threshold [σ]")
        plt.ylabel("Número de detecciones")
        plt.title("DAOFIND: detecciones vs. threshold")
        plt.legend()
        plt.tight_layout()
        plt.show()

        xgrid = np.linspace(
            thresholds.min(),
            thresholds.max(),
            300
        )

        plt.figure(figsize=(7, 4.5))
        plt.scatter(thresholds, lnN, label="Datos")
        plt.plot(
            xgrid,
            np.polyval([a, b, c], xgrid),
            label="Ajuste cuadrático"
        )
        plt.axvline(
            threshold_opt,
            linestyle="--",
            label=f"Threshold = {threshold_opt:.2f}"
        )
        plt.xlabel("Threshold [σ]")
        plt.ylabel(r"$\ln(N)$")
        plt.title(r"Ajuste de $\ln(N)$ vs. threshold")
        plt.legend()
        plt.tight_layout()
        plt.show()

    return resultado


# === FUNCIONES DEL NOTEBOOK, CELDA 116 ===
def cargar_daofind(archivo):
    columnas = [
        "x",
        "y",
        "mag_daofind",
        "sharpness",
        "sround",
        "ground",
        "id",
    ]

    return pd.read_csv(
        archivo,
        comment="#",
        sep=r"\s+",
        names=columnas,
        engine="python",
    )


def agregar_distancias_geometricas(
    imagen_path,
    catalogo,
    valor_invalido=VALOR_INVALIDO_DEFAULT,
):
    imagen = fits.getdata(imagen_path).astype(float)
    ny, nx = imagen.shape

    mask_invalidos = imagen == valor_invalido
    dist_invalidos = distance_transform_edt(~mask_invalidos)

    resultado = catalogo.copy()

    xpix = np.rint(resultado["x"]).astype(int).to_numpy() - 1
    ypix = np.rint(resultado["y"]).astype(int).to_numpy() - 1

    dentro = (
        (xpix >= 0) &
        (xpix < nx) &
        (ypix >= 0) &
        (ypix < ny)
    )

    d_invalidos = np.full(len(resultado), np.nan)

    d_invalidos[dentro] = dist_invalidos[
        ypix[dentro],
        xpix[dentro],
    ]

    x = resultado["x"].to_numpy()
    y = resultado["y"].to_numpy()

    d_borde = np.minimum.reduce([
        x - 1,
        nx - x,
        y - 1,
        ny - y,
    ])

    resultado["dist_invalidos"] = d_invalidos
    resultado["dist_borde"] = d_borde

    return resultado, imagen


def filtrar_distancia(catalogo, distancia):
    mask = (
        catalogo["dist_invalidos"].gt(distancia) &
        catalogo["dist_borde"].gt(distancia)
    )

    return (
        catalogo[mask].copy(),
        catalogo[~mask].copy(),
    )


# === FUNCIONES DEL NOTEBOOK, CELDA 118 ===
def cargar_phot(archivo):
    columnas = [
        "x",
        "y",
        "msky",
        "stdev",
        "nsky",
        "flux",
        "mag",
        "merr",
    ]

    return pd.read_csv(
        archivo,
        sep=r"\s+",
        names=columnas,
        na_values="INDEF",
    )


def seleccionar_calidad(df_phot, merr_max=MERR_MAX_DEFAULT):
    validas = df_phot[
        df_phot["mag"].notna() &
        df_phot["merr"].notna()
    ].copy()

    calidad = validas[
        validas["merr"] <= merr_max
    ].copy()

    return validas, calidad


def filtrar_vecinos(df_calidad, catalogo_daofind, d_crit):
    if df_calidad.empty:
        raise ValueError("No hay estrellas PHOT que pasen el corte de calidad.")

    targets = df_calidad[["x", "y"]].to_numpy()
    pool = catalogo_daofind[["x", "y"]].to_numpy()
    if len(pool) < 2:
        raise ValueError("Se requieren al menos dos detecciones DAOFIND para filtrar vecinos.")

    tree = cKDTree(pool)

    distancias, _ = tree.query(targets, k=2)

    d_match = distancias[:, 0]
    d_vecino = distancias[:, 1]

    resultado = df_calidad.copy()
    resultado["d_vecino"] = d_vecino

    mask = resultado["d_vecino"] > d_crit

    aisladas = resultado[mask].copy()
    rechazadas = resultado[~mask].copy()

    print(f"Máx. distancia match propio = {d_match.max():.6f}")
    print(f"Mediana distancia vecino    = {np.median(d_vecino):.4f} pix")
    print(f"Aisladas                    = {len(aisladas)}")
    print(f"Rechazadas por vecino       = {len(rechazadas)}")

    return aisladas, rechazadas


# === FUNCIONES DEL NOTEBOOK, CELDA 120 ===
def evaluar_anillo(
    imagen,
    x_iraf,
    y_iraf,
    r_in,
    r_out,
    valor_invalido=VALOR_INVALIDO_DEFAULT,
):
    ny, nx = imagen.shape

    x = float(x_iraf) - 1
    y = float(y_iraf) - 1

    if (
        x - r_out < 0 or
        x + r_out >= nx or
        y - r_out < 0 or
        y + r_out >= ny
    ):
        return None

    xmin = int(np.floor(x - r_out))
    xmax = int(np.ceil(x + r_out)) + 1
    ymin = int(np.floor(y - r_out))
    ymax = int(np.ceil(y + r_out)) + 1

    sub = imagen[ymin:ymax, xmin:xmax]

    yy, xx = np.indices(sub.shape)

    xxg = xx + xmin
    yyg = yy + ymin

    r = np.sqrt(
        (xxg - x) ** 2 +
        (yyg - y) ** 2
    )

    mask = (r >= r_in) & (r <= r_out)
    valores = sub[mask]

    if np.any(valores == valor_invalido):
        return None

    valores = valores[np.isfinite(valores)]

    if len(valores) == 0:
        return None

    q1 = np.percentile(valores, 25)
    q3 = np.percentile(valores, 75)
    iqr = q3 - q1
    limite = q3 + 1.5 * iqr

    n_sky = len(valores)
    n_out = int(np.sum(valores > limite))
    f_out = n_out / n_sky

    return {
        "Q1": q1,
        "Q3": q3,
        "IQR": iqr,
        "limite_alto": limite,
        "N_sky": n_sky,
        "N_out": n_out,
        "f_out": f_out,
    }


def analizar_anillo_fijo(
    imagen,
    df_aisladas,
    r_in,
    ancho=ANCHO_SKY_DEFAULT,
):
    r_out = r_in + ancho

    filas = []

    for _, estrella in df_aisladas.iterrows():
        res = evaluar_anillo(
            imagen,
            estrella["x"],
            estrella["y"],
            r_in,
            r_out,
        )

        if res is None:
            continue

        filas.append({
            "x": estrella["x"],
            "y": estrella["y"],
            **res,
        })

    df_anillos = pd.DataFrame(filas)
    if df_anillos.empty:
        raise ValueError("No se pudo medir ningún anillo de cielo fijo válido.")

    f = df_anillos["f_out"].to_numpy()

    q1 = np.percentile(f, 25)
    q3 = np.percentile(f, 75)
    iqr = q3 - q1
    f_lim = q3 + 1.5 * iqr

    df_limpio = df_anillos[
        df_anillos["f_out"] <= f_lim
    ].copy()

    df_contaminado = df_anillos[
        df_anillos["f_out"] > f_lim
    ].copy()

    print(f"Total analizadas = {len(df_anillos)}")
    print(f"f_lim             = {f_lim:.6f} ({100*f_lim:.2f}%)")
    print(f"Cielo limpio      = {len(df_limpio)}")
    print(f"Contaminadas      = {len(df_contaminado)}")

    return df_anillos, df_limpio, df_contaminado, float(f_lim)


def seleccionar_anillo_adaptativo(
    imagen,
    df_aisladas,
    f_lim,
    r_inicial,
    ancho=ANCHO_SKY_DEFAULT,
    paso=PASO_ANILLO_DEFAULT,
    n_anillos=N_ANILLOS_DEFAULT,
):
    radios = [
        (
            r_inicial + i * paso,
            r_inicial + i * paso + ancho,
        )
        for i in range(n_anillos)
    ]

    filas = []

    for _, estrella in df_aisladas.iterrows():
        x = estrella["x"]
        y = estrella["y"]

        solucion = None
        hubo_valido = False

        for r_in, r_out in radios:
            res = evaluar_anillo(
                imagen,
                x,
                y,
                r_in,
                r_out,
            )

            if res is None:
                continue

            hubo_valido = True

            if res["f_out"] <= f_lim:
                solucion = {
                    "x": x,
                    "y": y,
                    "r_in": r_in,
                    "r_out": r_out,
                    "f_out": res["f_out"],
                    "N_sky": res["N_sky"],
                    "N_out": res["N_out"],
                    "estado": "aceptado",
                }
                break

        if solucion is None:
            solucion = {
                "x": x,
                "y": y,
                "r_in": np.nan,
                "r_out": np.nan,
                "f_out": np.nan,
                "N_sky": np.nan,
                "N_out": np.nan,
                "estado": (
                    "contaminado_hasta_rmax"
                    if hubo_valido
                    else "sin_solucion_geometrica"
                ),
            }

        filas.append(solucion)

    return pd.DataFrame(filas)


# === FUNCIONES DEL NOTEBOOK, CELDA 122 ===
def crear_grupos_adaptativos(
    df_adapt,
    chip,
    filtro,
    base=BASE,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    prefijo = rutas["prefijo"]

    aceptadas = df_adapt[
        df_adapt["estado"] == "aceptado"
    ].copy()

    grupos = []

    for (r_in, r_out), grupo in aceptadas.groupby(
        ["r_in", "r_out"]
    ):
        tag = (
            f"r{r_in:.2f}_{r_out:.2f}"
            .replace(".", "_")
        )

        archivo = carpeta / f"{prefijo}_adapt_{tag}.coo"

        grupo[["x", "y"]].to_csv(
            archivo,
            sep=" ",
            index=False,
            header=False,
            float_format="%.3f",
        )

        grupos.append({
            "r_in": r_in,
            "r_out": r_out,
            "n": len(grupo),
            "coords": archivo,
            "tag": tag,
        })

    tabla = pd.DataFrame(grupos)

    print(tabla[["r_in", "r_out", "n"]])

    return tabla


def generar_script_phot_adaptativo(
    chip,
    filtro,
    tabla_grupos,
    fwhm,
    sigma_sky,
    base=BASE,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    imagen = rutas["imagen"]
    prefijo = rutas["prefijo"]

    hdr = fits.getheader(imagen)

    gain = hdr.get("GAIN", 2.75)
    rdnoise = hdr.get("RDNOISE", 6.4125)
    saturate = hdr.get("SATURATE", 39000)

    r_ap = 3 * fwhm

    lineas = [
        f'datapars.fwhmpsf = {fwhm:.6f}',
        f'datapars.sigma = {sigma_sky:.6f}',
        'datapars.datamin = INDEF',
        f'datapars.datamax = {saturate}',
        f'datapars.readnoise = {rdnoise}',
        f'datapars.epadu = {gain}',
        f'photpars.apertures = {r_ap:.6f}',
        f'fitskypars.dannulus = {ANCHO_SKY_DEFAULT}',
        'fitskypars.salgorithm = "mode"',
        'centerpars.calgorithm = "none"',
        '',
    ]

    for _, fila in tabla_grupos.iterrows():
        r_in = fila["r_in"]
        tag = fila["tag"]
        coords = fila["coords"]

        mag = carpeta / f"{prefijo}_adapt_{tag}.mag"
        dat = carpeta / f"{prefijo}_adapt_{tag}.dat"

        lineas += [
            f'fitskypars.annulus = {r_in:.6f}',
            (
                f'phot {imagen} coords={coords} '
                f'output={mag} interac- verify-'
            ),
            (
                f'txdump {mag} '
                '"XCENTER,YCENTER,MSKY,STDEV,NSKY,FLUX,MAG,MERR" '
                f'yes > {dat}'
            ),
            '',
        ]

    script = carpeta / f"{prefijo}_phot_adaptativo.cl"
    script.write_text(
        "\n".join(lineas) + "\n",
        encoding="utf-8"
    )

    print("Script creado:")
    print(script)
    print()
    print("Ejecutar dentro de IRAF:")
    print(f"cl < {script}")

    return script


# === FUNCIONES DEL NOTEBOOK, CELDA 132 ===
def _tag_num(x, nd=2):
    """Convierte 3.55 -> 3_55; 2.0 -> 2_0."""
    return f"{float(x):.{nd}f}".rstrip("0").rstrip(".").replace(".", "_")


def contar_catalogo_iraf(path):
    """Cuenta filas reales de un catálogo IRAF ignorando comentarios y blancos."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    n = 0
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for linea in f:
            s = linea.strip()
            if s and not s.startswith("#"):
                n += 1
    return n


def contar_pruebas_daofind(
    chip,
    filtro="V",
    thresholds=THRESHOLDS_DEFAULT,
    base=BASE,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    prefijo = rutas["prefijo"]

    filas = []

    for th in thresholds:
        # La función original usa str(th).replace(".", "_")
        tag_original = str(float(th)).replace(".", "_")
        archivo = carpeta / f"{prefijo}_thr{tag_original}.coo"

        # Compatibilidad con nombres sin .0
        if not archivo.exists():
            tag_alt = _tag_num(th)
            archivo_alt = carpeta / f"{prefijo}_thr{tag_alt}.coo"
            if archivo_alt.exists():
                archivo = archivo_alt

        filas.append({
            "threshold": float(th),
            "detecciones": contar_catalogo_iraf(archivo),
            "archivo": str(archivo),
        })

    tabla = pd.DataFrame(filas)
    print(tabla[["threshold", "detecciones"]].to_string(index=False))

    return tabla


# === FUNCIONES DEL NOTEBOOK, CELDA 134 ===
def etapa_threshold(
    chip,
    filtro,
    fwhm,
    sigma_sky,
    ejecutar_analisis=True,
):
    """
    1. Genera script de pruebas DAOFIND.
    2. Si los .coo ya existen, cuenta detecciones y calcula threshold.
    """
    script = generar_script_daofind(
        chip=chip,
        filtro=filtro,
        fwhm=fwhm,
        sigma_sky=sigma_sky,
    )

    rutas = rutas_chip(chip, filtro)
    prefijo = rutas["prefijo"]

    esperados = [
        rutas["carpeta"] / f"{prefijo}_thr{str(float(t)).replace('.', '_')}.coo"
        for t in THRESHOLDS_DEFAULT
    ]

    if not all(p.exists() for p in esperados):
        print()
        print("Todavía faltan archivos .coo.")
        print("Ejecuta primero en IRAF:")
        print(f"cl < {script}")
        return {"script": script, "tabla": None, "threshold": None}

    tabla = contar_pruebas_daofind(chip, filtro)

    if ejecutar_analisis:
        res = analizar_threshold(
            tabla["threshold"].to_numpy(),
            tabla["detecciones"].to_numpy(),
            graficar=True,
        )
        threshold = res["threshold_optimo"]
    else:
        res = None
        threshold = None

    return {
        "script": script,
        "tabla": tabla,
        "analisis": res,
        "threshold": threshold,
    }


# === FUNCIONES DEL NOTEBOOK, CELDA 136 ===
def generar_script_daofind_final(
    chip,
    filtro,
    fwhm,
    sigma_sky,
    threshold,
    base=BASE,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    imagen = rutas["imagen"]
    prefijo = rutas["prefijo"]

    hdr = fits.getheader(imagen)
    gain = hdr.get("GAIN", 2.75)
    rdnoise = hdr.get("RDNOISE", 6.4125)
    saturate = hdr.get("SATURATE", 39000)

    salida = carpeta / f"{prefijo}_final.coo"
    script = carpeta / f"{prefijo}_daofind_final.cl"

    lineas = [
        f"datapars.fwhmpsf = {fwhm:.6f}",
        f"datapars.sigma = {sigma_sky:.6f}",
        "datapars.datamin = INDEF",
        f"datapars.datamax = {saturate}",
        f"datapars.readnoise = {rdnoise}",
        f"datapars.epadu = {gain}",
        f"findpars.threshold = {threshold:.6f}",
        (
            f"daofind {imagen} "
            f"output={salida} verify-"
        ),
    ]

    script.write_text(
    "\n".join(lineas) + "\n",
    encoding="utf-8"
)

    print("Script DAOFIND final:")
    print(script)
    print()
    print("Ejecutar en IRAF:")
    print(f"cl < {script}")

    return script, salida


def preparar_phot_fijo(
    chip,
    filtro,
    fwhm,
    base=BASE,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    prefijo = rutas["prefijo"]
    imagen_path = rutas["imagen"]
    catalogo_path = carpeta / f"{prefijo}_final.coo"

    catalogo = cargar_daofind(catalogo_path)
    catalogo_geo, imagen = agregar_distancias_geometricas(
        imagen_path,
        catalogo,
    )

    r_ap = 3.0 * fwhm
    r_sky_in = 5.0 * fwhm
    r_sky_out = r_sky_in + ANCHO_SKY_DEFAULT

    catalogo_ap, rechazado_ap = filtrar_distancia(
        catalogo_geo,
        r_ap,
    )

    catalogo_phot, rechazado_sky = filtrar_distancia(
        catalogo_geo,
        r_sky_out,
    )

    xy = carpeta / f"{prefijo}_phot_ready.xy"
    catalogo_phot[["x", "y"]].to_csv(
        xy,
        sep=" ",
        index=False,
        header=False,
        float_format="%.3f",
    )

    diag = carpeta / f"{prefijo}_geometria.dat"
    catalogo_geo.to_csv(
        diag,
        sep=" ",
        index=False,
        float_format="%.5f",
    )

    print("DAOFIND inicial      :", len(catalogo_geo))
    print("Apertura válida      :", len(catalogo_ap))
    print("Anillo completo válido:", len(catalogo_phot))
    print()
    print("Archivo PHOT-ready:")
    print(xy)

    return {
        "catalogo": catalogo_geo,
        "catalogo_ap": catalogo_ap,
        "catalogo_phot": catalogo_phot,
        "imagen": imagen,
        "r_ap": r_ap,
        "r_sky_in": r_sky_in,
        "r_sky_out": r_sky_out,
        "xy": xy,
        "diag": diag,
    }


# === FUNCIONES DEL NOTEBOOK, CELDA 137 ===
def generar_script_phot_fijo(
    chip,
    filtro,
    fwhm,
    sigma_sky,
    base=BASE,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    imagen = rutas["imagen"]
    prefijo = rutas["prefijo"]

    hdr = fits.getheader(imagen)
    gain = hdr.get("GAIN", 2.75)
    rdnoise = hdr.get("RDNOISE", 6.4125)
    saturate = hdr.get("SATURATE", 39000)

    r_ap = 3.0 * fwhm
    r_sky_in = 5.0 * fwhm

    coords = carpeta / f"{prefijo}_phot_ready.xy"
    mag = carpeta / f"{prefijo}_fixed.mag"
    dat = carpeta / f"{prefijo}_fixed.dat"
    script = carpeta / f"{prefijo}_phot_fijo.cl"

    lineas = [
        'noao',
        'digiphot',
        'daophot',
        '',
        f"datapars.fwhmpsf = {fwhm:.6f}",
        f"datapars.sigma = {sigma_sky:.6f}",
        "datapars.datamin = INDEF",
        f"datapars.datamax = {saturate}",
        f"datapars.readnoise = {rdnoise}",
        f"datapars.epadu = {gain}",
        f"photpars.apertures = {r_ap:.6f}",
        f"fitskypars.annulus = {r_sky_in:.6f}",
        f"fitskypars.dannulus = {ANCHO_SKY_DEFAULT}",
        'fitskypars.salgorithm = "mode"',
        'centerpars.calgorithm = "none"',
        (
            f"phot {imagen} coords={coords} "
            f"output={mag} interac- verify-"
        ),
        (
            f'txdump {mag} '
            '"XCENTER,YCENTER,MSKY,STDEV,NSKY,FLUX,MAG,MERR" '
            f"yes > {dat}"
        ),
    ]

    script.write_text(
        "\n".join(lineas) + "\n",
        encoding="utf-8"
    )

    print("Script PHOT fijo:")
    print(script)
    print()
    print("Ejecutar en IRAF:")
    print(f"cl < {script}")

    return script, dat


# === FUNCIONES DEL NOTEBOOK, CELDA 139 ===
def analizar_post_phot(
    chip,
    filtro,
    fwhm,
    sigma_sky,
    base=BASE,
    merr_max=MERR_MAX_DEFAULT,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    prefijo = rutas["prefijo"]

    # Catálogo completo DAOFIND
    cat = cargar_daofind(carpeta / f"{prefijo}_final.coo")
    cat_geo, imagen = agregar_distancias_geometricas(
        rutas["imagen"],
        cat,
    )

    # PHOT fijo
    fixed_path = carpeta / f"{prefijo}_fixed.dat"
    df_fixed = cargar_phot(fixed_path)

    validas, calidad = seleccionar_calidad(
        df_fixed,
        merr_max=merr_max,
    )

    r_ap = 3.0 * fwhm
    r_sky_in = 5.0 * fwhm

    aisladas, rechazadas_vecino = filtrar_vecinos(
        calidad,
        cat_geo,
        d_crit=r_ap,
    )

    (
        df_anillos,
        cielo_limpio,
        cielo_contaminado,
        f_lim,
    ) = analizar_anillo_fijo(
        imagen,
        aisladas,
        r_in=r_sky_in,
        ancho=ANCHO_SKY_DEFAULT,
    )

    df_adapt = seleccionar_anillo_adaptativo(
        imagen,
        aisladas,
        f_lim=f_lim,
        r_inicial=r_sky_in,
        ancho=ANCHO_SKY_DEFAULT,
        paso=PASO_ANILLO_DEFAULT,
        n_anillos=N_ANILLOS_DEFAULT,
    )

    tabla_grupos = crear_grupos_adaptativos(
        df_adapt,
        chip,
        filtro,
        base=base,
    )

    # Guardar información esencial para recuperar el análisis
    calidad.to_csv(
        carpeta / f"{prefijo}_quality.dat",
        sep=" ",
        index=False,
        float_format="%.6f",
    )

    aisladas.to_csv(
        carpeta / f"{prefijo}_aisladas.dat",
        sep=" ",
        index=False,
        float_format="%.6f",
    )

    df_anillos.to_csv(
        carpeta / f"{prefijo}_anillo_fijo.dat",
        sep=" ",
        index=False,
        float_format="%.6f",
    )

    df_adapt.to_csv(
        carpeta / f"{prefijo}_anillo_adaptativo.dat",
        sep=" ",
        index=False,
        float_format="%.6f",
    )

    tabla_grupos.to_csv(
        carpeta / f"{prefijo}_grupos_adaptativos.csv",
        index=False,
    )

    script_adapt = generar_script_phot_adaptativo(
        chip,
        filtro,
        tabla_grupos,
        fwhm,
        sigma_sky,
        base=base,
    )

    n_limpias = len(cielo_limpio)
    n_contam = len(cielo_contaminado)
    n_aceptadas = (df_adapt["estado"] == "aceptado").sum()
    n_recuperadas = n_aceptadas - n_limpias

    resumen = {
        "chip": chip,
        "filtro": filtro,
        "fwhm": fwhm,
        "sigma_sky": sigma_sky,
        "N_phot_total": len(df_fixed),
        "N_phot_validas": len(validas),
        "N_MERR": len(calidad),
        "N_aisladas": len(aisladas),
        "f_lim": f_lim,
        "N_fijo_limpio": n_limpias,
        "N_fijo_contaminado": n_contam,
        "N_adapt_aceptadas": int(n_aceptadas),
        "N_recuperadas": int(n_recuperadas),
        "N_sin_solucion": int(len(df_adapt) - n_aceptadas),
        "tasa_recuperacion": (
            100.0 * n_recuperadas / n_contam
            if n_contam > 0 else np.nan
        ),
    }

    print()
    print("===== RESUMEN PRELIMINAR =====")
    for k, v in resumen.items():
        print(f"{k}: {v}")

    return {
        "df_fixed": df_fixed,
        "validas": validas,
        "calidad": calidad,
        "aisladas": aisladas,
        "rechazadas_vecino": rechazadas_vecino,
        "df_anillos": df_anillos,
        "cielo_limpio": cielo_limpio,
        "cielo_contaminado": cielo_contaminado,
        "df_adapt": df_adapt,
        "tabla_grupos": tabla_grupos,
        "script_adapt": script_adapt,
        "resumen": resumen,
    }


# === FUNCIONES DEL NOTEBOOK, CELDA 141 ===

def comparar_fijo_adaptativo(
    chip,
    filtro,
    r_sky_in,
    df_anillos=None,
    base=BASE,
):
    rutas = rutas_chip(chip, filtro, base)
    carpeta = rutas["carpeta"]
    prefijo = rutas["prefijo"]

    # ==========================================
    # 1. CARGAR FOTOMETRÍA ADAPTATIVA
    # ==========================================

    grupos_path = carpeta / f"{prefijo}_grupos_adaptativos.csv"
    tabla_grupos = pd.read_csv(grupos_path)

    cols = [
        "x", "y", "msky", "stdev",
        "nsky", "flux", "mag", "merr"
    ]

    grupos = []

    for _, fila in tabla_grupos.iterrows():
        tag = fila["tag"]
        archivo = carpeta / f"{prefijo}_adapt_{tag}.dat"

        temp = pd.read_csv(
            archivo,
            sep=r"\s+",
            names=cols,
            na_values=["INDEF"],
        )

        for col in cols:
            temp[col] = pd.to_numeric(
                temp[col], errors="coerce"
            )

        temp["r_in"] = float(fila["r_in"])
        temp["r_out"] = float(fila["r_out"])
        temp["tag"] = tag

        grupos.append(temp)

    if not grupos:
        raise ValueError("No hay grupos con fotometría adaptativa para comparar.")

    adapt = pd.concat(grupos, ignore_index=True)

    # ==========================================
    # 2. CARGAR FOTOMETRÍA FIJA
    # ==========================================

    fixed_path = carpeta / f"{prefijo}_fixed.dat"
    fixed = cargar_phot(fixed_path)
    if fixed.empty:
        raise ValueError("El catálogo PHOT fijo está vacío.")

    for col in cols:
        fixed[col] = pd.to_numeric(
            fixed[col], errors="coerce"
        )

    # ==========================================
    # 3. ASOCIAR FUENTES POR COORDENADAS
    # ==========================================

    if adapt[["x", "y"]].isna().any().any():
        raise ValueError(
            "Hay coordenadas inválidas en PHOT adaptativo."
        )

    if fixed[["x", "y"]].isna().any().any():
        raise ValueError(
            "Hay coordenadas inválidas en PHOT fijo."
        )

    tree = cKDTree(
        fixed[["x", "y"]].to_numpy()
    )

    dmatch, idx = tree.query(
        adapt[["x", "y"]].to_numpy(),
        k=1
    )

    max_dist = float(np.max(dmatch))

    print("===== ASOCIACIÓN DE FUENTES =====")
    print("Matches:", len(adapt))
    print(f"Máx. distancia match: {max_dist:.6f} pix")

    # Evitar asociaciones de estrellas distintas
    if max_dist > 0.1:
        raise ValueError(
            "Existen asociaciones con distancia > 0.1 pix. "
            "Revisar coordenadas antes de comparar."
        )

    comp = adapt.rename(
        columns={
            col: f"{col}_adapt"
            for col in ["msky", "stdev", "flux", "mag", "merr"]
        }
    ).copy()

    comp["dist_match"] = dmatch

    for col in ["msky", "stdev", "flux", "mag", "merr"]:
        comp[f"{col}_fixed"] = (
            fixed.iloc[idx][col].to_numpy()
        )

    # ==========================================
    # 4. VALIDACIÓN DEL GRUPO DE CONTROL
    # ==========================================

    control = comp[
        np.isclose(comp["r_in"], r_sky_in)
    ].copy()

    print()
    print("===== GRUPO DE CONTROL =====")
    print("Estrellas de control:", len(control))

    for col in ["mag", "merr", "msky", "stdev", "flux"]:
        delta = (
            control[f"{col}_adapt"]
            - control[f"{col}_fixed"]
        )

        control[f"delta_{col}"] = delta

        valores = delta.abs().dropna()

        max_delta = (
            valores.max()
            if len(valores)
            else np.nan
        )

        print(f"Máx |Δ{col.upper()}|: {max_delta}")

    # ==========================================
    # 5. ESTRELLAS CON ANILLO DESPLAZADO
    # ==========================================

    movidas = comp[
        comp["r_in"] > r_sky_in + 1e-8
    ].copy()

    n_movidas_total = len(movidas)

    required = [
        "mag_fixed", "mag_adapt",
        "merr_fixed", "merr_adapt",
        "msky_fixed", "msky_adapt",
        "stdev_fixed", "stdev_adapt",
        "flux_fixed", "flux_adapt",
    ]

    movidas = movidas.dropna(
        subset=required
    ).copy()

    # Evitar divisiones por cero o valores no físicos
    movidas = movidas[
        (movidas["merr_fixed"] > 0) &
        (movidas["merr_adapt"] > 0) &
        (movidas["flux_fixed"] > 0) &
        (movidas["flux_adapt"] > 0)
    ].copy()

    print()
    print("===== ANILLOS DESPLAZADOS =====")
    print("Total desplazadas:", n_movidas_total)
    print("Fotometría válida en ambas:", len(movidas))

    # ==========================================
    # 6. DIFERENCIAS FOTOMÉTRICAS
    # ==========================================

    for col in ["msky", "stdev", "mag", "merr"]:
        movidas[f"delta_{col}"] = (
            movidas[f"{col}_adapt"]
            - movidas[f"{col}_fixed"]
        )

    movidas["flux_ratio"] = (
        movidas["flux_adapt"]
        / movidas["flux_fixed"]
    )

    movidas["snr_fixed"] = (
        1.0857 / movidas["merr_fixed"]
    )

    movidas["snr_adapt"] = (
        1.0857 / movidas["merr_adapt"]
    )

    # ==========================================
    # 7. RESUMEN DE COMPARACIÓN
    # ==========================================

    def mediana(col):
        if len(movidas) == 0:
            return np.nan
        return float(np.nanmedian(movidas[col]))

    resumen = {
        "N_control": len(control),
        "N_movidas_total": n_movidas_total,
        "N_movidas_validas": len(movidas),

        "med_delta_msky": mediana("delta_msky"),
        "med_delta_stdev": mediana("delta_stdev"),
        "med_delta_mag": mediana("delta_mag"),

        "med_abs_delta_mag": (
            float(np.nanmedian(np.abs(movidas["delta_mag"])))
            if len(movidas) else np.nan
        ),

        "med_merr_fixed": mediana("merr_fixed"),
        "med_merr_adapt": mediana("merr_adapt"),

        "pct_merr_mejora": (
            100 * np.mean(
                movidas["merr_adapt"]
                < movidas["merr_fixed"]
            )
            if len(movidas) else np.nan
        ),

        "med_snr_fixed": mediana("snr_fixed"),
        "med_snr_adapt": mediana("snr_adapt"),
        "med_flux_ratio": mediana("flux_ratio"),
    }

    # ==========================================
    # 8. ANALIZAR CONTAMINACIÓN DEL ANILLO
    # ==========================================

    if df_anillos is None:
        path_anillos = (
            carpeta / f"{prefijo}_anillo_fijo.dat"
        )
        df_anillos = pd.read_csv(
            path_anillos,
            sep=r"\s+"
        )

    path_adapt = (
        carpeta / f"{prefijo}_anillo_adaptativo.dat"
    )

    df_sel = pd.read_csv(
        path_adapt,
        sep=r"\s+",
        na_values=["INDEF"]
    )

    # CORRECCIÓN PRINCIPAL DEL ERROR TypeError
    for col in ["x", "y", "r_in", "f_out"]:
        df_sel[col] = pd.to_numeric(
            df_sel[col],
            errors="coerce"
        )

    aceptadas = df_sel["estado"] == "aceptado"

    assert df_sel.loc[aceptadas, "r_in"].notna().all(), (
        "Existen estrellas aceptadas sin r_in válido."
    )

    assert df_sel.loc[aceptadas, "f_out"].notna().all(), (
        "Existen estrellas aceptadas sin f_out válido."
    )

    print()
    print("===== SELECCIÓN ADAPTATIVA =====")
    print("Estrellas aceptadas:", int(aceptadas.sum()))

    # Fuentes que realmente cambiaron de anillo
    rec = df_sel[
        aceptadas &
        (df_sel["r_in"] > r_sky_in + 1e-8)
    ][["x", "y", "f_out"]].copy()

    rec = rec.rename(
        columns={"f_out": "f_out_final"}
    )

    # Contaminación del anillo original
    ini = df_anillos[
        ["x", "y", "f_out"]
    ].copy()

    ini = ini.rename(
        columns={"f_out": "f_out_inicial"}
    )

    for col in ["x", "y", "f_out_inicial"]:
        ini[col] = pd.to_numeric(
            ini[col],
            errors="coerce"
        )

    # Redondear coordenadas para evitar errores
    # por diferencias de precisión al guardar archivos
    for tabla in [rec, ini]:
        tabla["x_key"] = tabla["x"].round(3)
        tabla["y_key"] = tabla["y"].round(3)

    claves = ["x_key", "y_key"]

    if rec.duplicated(claves).any():
        raise ValueError(
            "Coordenadas duplicadas en anillos desplazados."
        )

    if ini.duplicated(claves).any():
        raise ValueError(
            "Coordenadas duplicadas en anillos iniciales."
        )

    fout = rec.merge(
        ini[claves + ["f_out_inicial"]],
        on=claves,
        how="left",
        validate="one_to_one"
    )

    if fout["f_out_inicial"].isna().any():
        raise ValueError(
            "No se pudo recuperar f_out inicial "
            "para algunas estrellas desplazadas."
        )

    fout["delta_f_out"] = (
        fout["f_out_final"]
        - fout["f_out_inicial"]
    )

    # ==========================================
    # 9. ESTADÍSTICAS DE CONTAMINACIÓN
    # ==========================================

    resumen["med_fout_inicial"] = (
        float(np.nanmedian(fout["f_out_inicial"]))
        if len(fout) else np.nan
    )

    resumen["med_fout_final"] = (
        float(np.nanmedian(fout["f_out_final"]))
        if len(fout) else np.nan
    )

    resumen["med_delta_fout"] = (
        float(np.nanmedian(fout["delta_f_out"]))
        if len(fout) else np.nan
    )

    resumen["pct_fout_reduce"] = (
        100 * np.mean(fout["delta_f_out"] < 0)
        if len(fout) else np.nan
    )

    # ==========================================
    # 10. RESUMEN FINAL
    # ==========================================

    print()
    print("===== COMPARACIÓN FOTOMÉTRICA =====")

    for k, v in resumen.items():
        print(f"{k}: {v}")

    return {
        "comparacion": comp,
        "control": control,
        "movidas": movidas,
        "fout": fout,
        "resumen": resumen,
    }


NC, NR, NBOX = 5, 3, 15

def _cfgpaths(chip, filtro, base):
    r = rutas_chip(chip, filtro, base)
    if not r['imagen'].exists():
        raise FileNotFoundError(f'FITS inexistente: {r["imagen"]}')
    return r


def _cat(path):
    c = cargar_daofind(path)
    for k in ('x','y','sharpness','sround','ground'):
        c[k] = pd.to_numeric(c[k], errors='coerce')
    return c.dropna(subset=['x','y'])


def _cell(x,y,shape):
    ny,nx=shape
    return min(NR-1,int(NR*y/ny)),min(NC-1,int(NC*x/nx))


def _skybox(image,x,y):
    """Medición de una caja 15x15, centro x/y cero-based."""
    ny,nx=image.shape;x=int(x);y=int(y); h=NBOX//2
    if min(x-h,y-h,nx-1-x-h,ny-1-y-h)<0:return None
    a=image[y-h:y+h+1,x-h:x+h+1]
    if a.shape!=(NBOX,NBOX) or not np.isfinite(a).all() or (a==50000).any():return None
    med=float(np.median(a));sigma=float(a.std(ddof=1))
    mad=float(1.4826*np.median(np.abs(a-med)))
    if sigma<=0 or mad<=0:return None
    q1,q3=np.percentile(a,[25,75]); high=q3+1.5*(q3-q1)
    # Sharpness inverso: fondo con escasa estructura positiva compacta.
    # No es el SHARP numérico de DAOFIND: éste solo existe en máximos detectados.
    sharp_proxy=float(np.max(a-gaussian_filter(a,1.1))/mad)
    return dict(x=x+1,y=y+1,sigma=sigma,mediana=med,
                f_out=float(np.mean(a>high)),sharp_proxy=sharp_proxy,
                pico_sigma=float((a.max()-med)/mad))


def _init_sky(image,seed):
    ny,nx=image.shape;rng=np.random.default_rng(seed);vals=[]
    for i in range(NR):
        for j in range(NC):
            group=[]
            for _ in range(80):
                x=int(rng.integers(j*nx//NC+16,(j+1)*nx//NC-16))
                y=int(rng.integers(i*ny//NR+16,(i+1)*ny//NR-16))
                d=_skybox(image,x,y)
                if d and d['pico_sigma']<5.5 and d['f_out']<.05:group.append(d['sigma'])
            if group:vals.append(np.percentile(group,25))
    if len(vals)<10:raise RuntimeError('No hay suficientes cajas para sigma inicial; revisar FITS.')
    return float(np.median(vals))


def _fit_fwhm(image,x,y,seed_fwhm,sigma,sat):
    """Ajuste gaussiano circular 2D; no idéntico al DIRECT de imexamine."""
    cx,cy=float(x)-1,float(y)-1;h=int(math.ceil(seed_fwhm*2.6))
    ix,iy=round(cx),round(cy);ny,nx=image.shape
    if min(ix-h,iy-h,nx-1-ix-h,ny-1-iy-h)<0:return None
    a=image[iy-h:iy+h+1,ix-h:ix+h+1]
    if not np.isfinite(a).all() or (a==50000).any() or a.max()>.85*sat:return None
    yy,xx=np.indices(a.shape,dtype=float);xx+=ix-h;yy+=iy-h
    r=np.hypot(xx-cx,yy-cy); ann=(r>seed_fwhm*1.8)&(r<seed_fwhm*2.5)
    if ann.sum()<20:return None
    bg=float(np.median(a[ann]));peak=float(a.max()-bg)
    if peak<max(12*sigma,25):return None
    focus=(r<1.5*seed_fwhm);X,Y,Z=xx[focus],yy[focus],a[focus]
    init=seed_fwhm/2.35482
    def res(p):
        amp,x0,y0,sig,b=p
        return (b+amp*np.exp(-.5*((X-x0)**2+(Y-y0)**2)/sig**2)-Z)/max(sigma,1)
    try:
        f=least_squares(res,[peak,cx,cy,init,bg],
            bounds=([peak*.1,cx-1.5,cy-1.5,init*.65,bg-4*sigma],
                    [peak*2.1,cx+1.5,cy+1.5,init*1.65,bg+4*sigma]),max_nfev=75)
    except (ValueError,FloatingPointError):return None
    if not f.success or np.hypot(f.x[1]-cx,f.x[2]-cy)>1.0:return None
    width=float(2.35482*f.x[3]);err=float(np.sqrt(np.mean(res(f.x)**2)))
    if not (.7*seed_fwhm<width<1.6*seed_fwhm) or err>6:return None
    return {'fwhm':width,'pico':float(a.max()),'rms_fit':err}


def _stars15(image,cat,guess,sigma,sat):
    ny,nx=image.shape
    good=cat.dropna(subset=['sharpness','sround','ground']).copy()
    good=good[(good.sharpness.between(.2,1.0)) &
              (good.sround.abs()<.5)&(good.ground.abs()<.5)].copy()
    if len(good)<NBOX:raise RuntimeError('No hay 15 detecciones con SHARP de tipo estelar.')
    tree=cKDTree(cat[['x','y']].to_numpy())
    bad=distance_transform_edt(np.isfinite(image)&(image!=50000))
    cells={(i,j):[] for i in range(NR) for j in range(NC)}
    for _,s in good.iterrows():
        x,y=float(s.x),float(s.y);ix,iy=round(x-1),round(y-1)
        if not (0<=ix<nx and 0<=iy<ny):continue
        if min(x-1,nx-x,y-1,ny-y)<4*guess or bad[iy,ix]<4*guess:continue
        if image[iy,ix]>.85*sat:continue
        d,_=tree.query([x,y],k=2)
        if d[0]>.05 or d[1]<5*guess:continue
        cell=_cell(x-1,y-1,image.shape)
        score=2*abs(s.sharpness-.6)+abs(s.sround)+abs(s.ground)
        cells[cell].append((score,s))
    found=[];missing=[]
    for c,values in cells.items():
        chosen=None
        for _,s in sorted(values,key=lambda x:x[0])[:35]:
            fit=_fit_fwhm(image,s.x,s.y,guess,sigma,sat)
            if fit:
                chosen=dict(x=float(s.x),y=float(s.y),sharpness=float(s.sharpness),
                            sround=float(s.sround),ground=float(s.ground),
                            fila=c[0]+1,columna=c[1]+1,**fit)
                break
        if chosen is None:missing.append(c)
        else:found.append(chosen)
    if missing:raise RuntimeError(f'Faltan estrellas aptas en celdas 5x3 {missing}. Inspeccionar o modificar FWHM inicial.')
    result=pd.DataFrame(found)
    if not 1<float(result.fwhm.median())<15:raise RuntimeError('FWHM automático no plausible.')
    return result


def _background15(image,cat,fwhm,sig0,seed):
    ny,nx=image.shape;rng=np.random.default_rng(seed)
    tree=cKDTree(cat[['x','y']].to_numpy());chosen=[];missing=[]
    for i in range(NR):
        for j in range(NC):
            trials=[]
            for _ in range(300):
                x=int(rng.integers(j*nx//NC+20,(j+1)*nx//NC-20))
                y=int(rng.integers(i*ny//NR+20,(i+1)*ny//NR-20))
                if tree.query([x+1,y+1])[0]<4*fwhm+11:continue
                d=_skybox(image,x,y)
                if d is None:continue
                if d['pico_sigma']>4.5 or d['sharp_proxy']>3.6 or d['f_out']>.06:continue
                if not (.45*sig0<d['sigma']<2*sig0):continue
                d.update(fila=i+1,columna=j+1)
                d['score']=d['sharp_proxy']+15*d['f_out']+.15*abs(d['sigma']/sig0-1)
                trials.append(d)
            if trials:chosen.append(min(trials,key=lambda t:t['score']))
            else:missing.append((i+1,j+1))
    if missing:raise RuntimeError(f'Faltan cajas de cielo en celdas {missing}. Revisar estructuras y fuentes.')
    d=pd.DataFrame(chosen);q1,q3=np.percentile(d.sigma,[25,75]);lim=q3+1.5*(q3-q1)
    d['aceptada_tukey']=d.sigma<=lim;d['limite_tukey']=lim
    if d.aceptada_tukey.sum()<12:raise RuntimeError('Menos de 12 de 15 regiones del cielo pasan Tukey.')
    return d


def _write_json(p,obj):
    Path(p).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')


def _save_plot(image,stars,sky,path,title):
    ny,nx=image.shape;scale=max(1,math.ceil(max(nx,ny)/1200))
    show=image[::scale,::scale]
    vals=show[np.isfinite(show)&(show!=50000)]
    lo,hi=np.percentile(vals,[5,99.5])
    fig,ax=plt.subplots(figsize=(12,7))
    ax.imshow(np.clip(show,lo,hi),origin='lower',cmap='gray',vmin=lo,vmax=hi,
              extent=[0,nx,0,ny],aspect='auto')
    ax.scatter(stars.x-1,stars.y-1,edgecolors='lime',facecolors='none',s=90,label='15 estrellas SHARP')
    ax.scatter(sky.x-1,sky.y-1,edgecolors='cyan',facecolors='none',s=95,marker='s',label='15 cajas fondo')
    for x in range(1,NC):ax.axvline(nx*x/NC,alpha=.25,lw=.5,color='white')
    for y in range(1,NR):ax.axhline(ny*y/NR,alpha=.25,lw=.5,color='white')
    ax.set(title=title,xlabel='x [pix]',ylabel='y [pix]');ax.legend()
    fig.tight_layout();fig.savefig(path,dpi=135);plt.close(fig)


def _prepare(r,args):
    im=fits.getdata(r['imagen']).astype(float)
    cat=_cat(r['carpeta']/f'{r["prefijo"]}_bootstrap.coo')
    if len(cat)<100:raise RuntimeError('Catálogo preliminar de DAOFIND muy pequeño.')
    sig0=args.seed_sigma if args.seed_sigma is not None else _init_sky(im,args.seed)
    header=fits.getheader(r['imagen'])
    sat=float(header.get('SATURATE',39000))
    stars=_stars15(im,cat,args.seed_fwhm,sig0,sat)
    # Control robusto adicional: FWHM extremos se documentan y no pesan en la mediana.
    q1, q3 = np.percentile(stars['fwhm'], [25, 75])
    width = q3-q1
    stars['aceptada_fwhm'] = stars['fwhm'].between(q1-1.5*width, q3+1.5*width)
    if int(stars['aceptada_fwhm'].sum()) < 12:
        raise RuntimeError('Menos de 12 de 15 ajustes de FWHM pasan el control de Tukey.')
    fwhm=round(float(stars.loc[stars.aceptada_fwhm, 'fwhm'].median()),2)
    sky=_background15(im,cat,fwhm,sig0,args.seed+71)
    sigma=round(float(sky.loc[sky.aceptada_tukey,'sigma'].median()),2)
    pref=r['prefijo'];dst=r['carpeta']
    s_csv=dst/f'{pref}_auto_estrellas.csv';sky_csv=dst/f'{pref}_auto_cielo.csv'
    stars.to_csv(s_csv,index=False,float_format='%.6f');sky.to_csv(sky_csv,index=False,float_format='%.6f')
    pic=dst/f'{pref}_auto_control.png'
    _save_plot(im,stars,sky,pic,f'{pref} | FWHM={fwhm:.2f} pix | σ cielo={sigma:.2f} ADU')
    info=dict(fwhm=fwhm,sigma_sky=sigma,sigma_inicial=sig0,n_estrellas=len(stars),
              n_cajas=len(sky),n_estrellas_fwhm_tukey=int(stars.aceptada_fwhm.sum()),
              n_cajas_tukey=int(sky.aceptada_tukey.sum()),
              stars_csv=str(s_csv),sky_csv=str(sky_csv),grafica=str(pic),
              metodo_fwhm='Ajuste gaussiano 2D a 15 fuentes con SHARP estelar de IRAF',
              metodo_cielo='15 cajas 15x15 sin fuentes y con proxy de sharpness bajo')
    _write_json(dst/f'{pref}_auto_config.json',info)
    print('15 estrellas distribuidas:',len(stars),'| 15 cajas cielo:',len(sky))
    print('FWHM:',fwhm,'pix | sigma del cielo:',sigma,'ADU')
    print('CONTROL VISUAL IMPORTANTE:',pic)
    return info


def _params(r):
    path=r['carpeta']/f'{r["prefijo"]}_auto_config.json'
    if not path.exists():raise RuntimeError(f'No hay FWHM/sigma automáticos: {path}. Ejecuta --stage select.')
    return json.loads(path.read_text())


def _iraf_params(r,fw,sigma):
    hdr=fits.getheader(r['imagen'])
    return [
      f'datapars.fwhmpsf = {fw:.6f}',f'datapars.sigma = {sigma:.6f}',
      'datapars.datamin = INDEF',
      f'datapars.datamax = {float(hdr.get("SATURATE",39000))}',
      f'datapars.readnoise = {float(hdr.get("RDNOISE",6.4125))}',
      f'datapars.epadu = {float(hdr.get("GAIN",2.75))}',
      ''
    ]


def _write_cl(path,lines):
    # Precarga paquetes y SIEMPRE termina con salto de línea.
    lines=['noao','digiphot','daophot','']+list(lines)
    path.write_text('\n'.join(lines).rstrip()+'\n',encoding='utf-8')
    return path


def _iraf_exec(script,outputs,args):
    outputs=[Path(p) for p in outputs]
    exist=[p.exists() for p in outputs]
    if all(exist) and not args.force:
        print(f'Ya existen los {len(outputs)} productos; no se recalculan.')
        return True
    if any(exist) and not args.force:
        raise RuntimeError('Salidas IRAF parciales existentes. Revisar y repetir con --force para hacer copia .bak y recalcular.')
    if args.mode=='manual':
        print(f'\nEJECUTAR EN IRAF: cl < {script}\nLuego invoca siguiente --stage.\n')
        return False
    if not shutil.which(args.iraf_exe):
        raise RuntimeError(f'No encuentro {args.iraf_exe}. Usa --iraf-exe ecl o --mode manual.')
    for p in outputs:
        if p.exists():
            shutil.copy2(p,p.with_suffix(p.suffix+'.bak'));p.unlink()
    log=script.with_suffix('.log')
    print('Ejecutando:',args.iraf_exe,'-f',script)
    with log.open('w') as f:
        proc=subprocess.run([args.iraf_exe,'-f',str(script)],cwd=args.base,
            stdout=f,stderr=subprocess.STDOUT,input='logout\n',text=True,
            timeout=args.timeout)
    if proc.returncode or not all(p.exists() for p in outputs):
        tail='\n'.join(log.read_text(errors='replace').splitlines()[-25:])
        raise RuntimeError(f'IRAF falló. Revisar {log}\n{tail}')
    txt=log.read_text(errors='replace')
    if 'INTERNAL ERROR' in txt or 'ERROR:' in txt:
        raise RuntimeError(f'IRAF reportó error; revisar log {log}')
    print('IRAF OK; log:',log)
    return True


def _need(*paths):
    for p in paths:
        if not Path(p).exists():raise RuntimeError(f'Falta archivo de fase anterior: {p}')


def _boot(r,args):
    im=fits.getdata(r['imagen']).astype(float)
    sig=args.seed_sigma if args.seed_sigma is not None else _init_sky(im,args.seed)
    cat=r['carpeta']/f'{r["prefijo"]}_bootstrap.coo'
    script=r['carpeta']/f'{r["prefijo"]}_bootstrap.cl'
    lines=_iraf_params(r,args.seed_fwhm,sig)
    lines+=['findpars.threshold = 5.0',f'daofind {r["imagen"]} output={cat} verify-']
    _write_cl(script,lines)
    return _iraf_exec(script,[cat],args)


def _thresholds(r,args):
    p=_params(r);ts=[2.,2.5,3.,3.5,4.,5.,6.]
    paths=[r['carpeta']/f'{r["prefijo"]}_thr{str(t).replace(".","_")}.coo' for t in ts]
    lines=_iraf_params(r,p['fwhm'],p['sigma_sky'])
    for t,o in zip(ts,paths):
        lines.extend([f'findpars.threshold = {t:.2f}',f'daofind {r["imagen"]} output={o} verify-',''])
    script=r['carpeta']/f'{r["prefijo"]}_daofind.cl'
    _write_cl(script,lines)
    completed=_iraf_exec(script,paths,args)
    if not completed:return False
    ns=[len(_cat(o)) for o in paths]
    if min(ns)<=0:raise RuntimeError('Un threshold tiene 0 detecciones.')
    fit=analizar_threshold(ts,ns,graficar=False)
    raw=fit['threshold_optimo'];th=round(raw,2)
    if not 2<=th<=6:raise RuntimeError(f'Threshold fuera del rango probado [2,6]: {th}. Revisión necesaria.')
    result={'threshold':th,'threshold_ajustado':raw,'r2':fit['r2'],'conteos':dict(zip(map(str,ts),ns))}
    _write_json(r['carpeta']/f'{r["prefijo"]}_threshold.json',result)
    pd.DataFrame({'threshold':ts,'N':ns}).to_csv(r['carpeta']/f'{r["prefijo"]}_thresholds.csv',index=False)
    fig,ax=plt.subplots(figsize=(7,4));ax.plot(ts,ns,'o-');ax.axvline(th,ls='--',label=f'th = {th:.2f}')
    ax.set(xlabel='Threshold [σ]',ylabel='N detectadas',title=r['prefijo']);ax.legend()
    fig.tight_layout();fig.savefig(r['carpeta']/f'{r["prefijo"]}_threshold.png');plt.close(fig)
    return True


def _final(r,args):
    p=_params(r);cfg=r['carpeta']/f'{r["prefijo"]}_threshold.json';_need(cfg)
    th=json.loads(cfg.read_text())['threshold']
    out=r['carpeta']/f'{r["prefijo"]}_final.coo'
    s=r['carpeta']/f'{r["prefijo"]}_daofind_final.cl'
    lines=_iraf_params(r,p['fwhm'],p['sigma_sky'])+[
        f'findpars.threshold = {th:.6f}',f'daofind {r["imagen"]} output={out} verify-']
    _write_cl(s,lines)
    return _iraf_exec(s,[out],args)


def _fixed(r,args):
    p=_params(r);_need(r['carpeta']/f'{r["prefijo"]}_final.coo')
    # Conserva exactamente los criterios y archivos del notebook.
    geo=preparar_phot_fijo(args.chip,args.filtro,p['fwhm'],base=args.base)
    script,dat=generar_script_phot_fijo(args.chip,args.filtro,p['fwhm'],p['sigma_sky'],base=args.base)
    # No duplicar encabezado si ya existía en la función original.
    content=script.read_text()
    if not content.startswith('noao\n'):content='noao\ndigiphot\ndaophot\n\n'+content
    script.write_text(content.rstrip()+'\n')
    mag=r['carpeta']/f'{r["prefijo"]}_fixed.mag'
    complete=_iraf_exec(script,[dat,mag],args)
    if complete:
        n0=contar_catalogo_iraf(geo['xy']);n1=contar_catalogo_iraf(dat)
        if n0!=n1:raise RuntimeError(f'PHOT fijo incompleto {n0} coordenadas vs {n1} mediciones.')
        print('PHOT fijo:',n0,'coordenadas =',n1,'mediciones.')
    return complete


def _adaptive(r,args):
    p=_params(r);_need(r['carpeta']/f'{r["prefijo"]}_fixed.dat')
    info=analizar_post_phot(args.chip,args.filtro,p['fwhm'],p['sigma_sky'],base=args.base)
    script=info['script_adapt'];content=script.read_text()
    if not content.startswith('noao\n'):content='noao\ndigiphot\ndaophot\n\n'+content
    script.write_text(content.rstrip()+'\n')
    outs=[]
    for _,row in info['tabla_grupos'].iterrows():
        for ext in ['mag','dat']:
            outs.append(r['carpeta']/f'{r["prefijo"]}_adapt_{row["tag"]}.{ext}')
    if not outs:raise RuntimeError('No se encontraron grupos de cielo aceptados.')
    complete=_iraf_exec(script,outs,args)
    if complete:
        for _,row in info['tabla_grupos'].iterrows():
            dat=r['carpeta']/f'{r["prefijo"]}_adapt_{row["tag"]}.dat'
            if contar_catalogo_iraf(dat)!=int(row['n']):
                raise RuntimeError(f'PHOT adaptativo incompleto: {dat}')
    return complete


def _compare(r,args):
    p=_params(r)
    groups=r['carpeta']/f'{r["prefijo"]}_grupos_adaptativos.csv';_need(groups)
    g=pd.read_csv(groups)
    _need(*[r['carpeta']/f'{r["prefijo"]}_adapt_{row["tag"]}.dat' for _,row in g.iterrows()])
    compdata=comparar_fijo_adaptativo(args.chip,args.filtro,5*p['fwhm'],base=args.base)
    comp=compdata['comparacion'];control=compdata['control'];moved=compdata['movidas'];fout=compdata['fout']
    if len(comp)!=g['n'].sum():raise RuntimeError('Falta alguna asociación en comparación final.')
    if len(comp) and comp['dist_match'].max()>0.001:
        raise RuntimeError('Asociación espacial no exacta: revisar coordenadas de PHOT.')
    if len(control)!=g.loc[np.isclose(g.r_in,5*p['fwhm']),'n'].sum():
        raise RuntimeError('Grupo control de tamaño incorrecto.')
    if control.empty:
        raise RuntimeError('No hay estrellas de control para validar PHOT adaptativo.')
    for var in ['mag','merr','msky','stdev','flux']:
        delta=(control[f'{var}_adapt']-control[f'{var}_fixed']).abs().dropna()
        if len(delta) != len(control) or delta.max()>1e-6:
            raise RuntimeError(f'El control cambió para {var}, no aceptar comparación adaptativa.')
    comp.to_csv(r['carpeta']/f'{r["prefijo"]}_comparacion.csv',index=False)
    moved.to_csv(r['carpeta']/f'{r["prefijo"]}_movidas.csv',index=False)
    fout.to_csv(r['carpeta']/f'{r["prefijo"]}_contaminacion.csv',index=False)
    result={'chip':args.chip,'filtro':args.filtro,**p,**compdata['resumen']}
    _write_json(r['carpeta']/f'{r["prefijo"]}_resumen_final.json',result)
    if len(moved):
        fig,ax=plt.subplots(figsize=(7,4))
        ax.hist(moved.delta_mag,bins=16);ax.axvline(0,ls='--')
        ax.set(xlabel='MAG adaptativa − fija',ylabel='N',title=r['prefijo']+' Cambio de magnitud')
        fig.tight_layout();fig.savefig(r['carpeta']/f'{r["prefijo"]}_delta_mag.png',dpi=130);plt.close(fig)
    if len(fout):
        fig,ax=plt.subplots(figsize=(6,5));ax.scatter(fout.f_out_inicial*100,fout.f_out_final*100)
        ax.set(xlabel='f_out inicial [%]',ylabel='f_out final [%]',title=r['prefijo']+' Contaminación')
        fig.tight_layout();fig.savefig(r['carpeta']/f'{r["prefijo"]}_delta_fout.png',dpi=130);plt.close(fig)
    print('Comparación validada; control:',len(control),'| fuentes desplazadas válidas:',len(moved))
    return True


def _stage(stage,r,args):
    return {
        'bootstrap':lambda:_boot(r,args),
        'select':lambda:_prepare(r,args),
        'thresholds':lambda:_thresholds(r,args),
        'final':lambda:_final(r,args),
        'fixed':lambda:_fixed(r,args),
        'adaptive':lambda:_adaptive(r,args),
        'compare':lambda:_compare(r,args),
    }[stage]()


def main():
    parser=argparse.ArgumentParser(
        description='Sextans A: 15 estrellas seleccionadas por SHARP y 15 regiones de cielo distribuidas; DAOFIND + PHOT + anillo adaptativo.')
    parser.add_argument('--base',type=Path,default=BASE,help='Directorio dataSextansA del proyecto')
    parser.add_argument('--chip',type=int,required=True,choices=range(1,9))
    parser.add_argument('--filtro',choices=['V','B','I'],default='V')
    parser.add_argument('--stage',choices=['all','bootstrap','select','thresholds','final','fixed','adaptive','compare'],default='all')
    parser.add_argument('--mode',choices=['auto','manual'],default='auto')
    parser.add_argument('--iraf-exe',default='irafcl',help='irafcl, ecl o cl en macOS')
    parser.add_argument('--seed-fwhm',type=float,default=4.5,help='FWHM inicial aproximado; no reemplaza ajuste')
    parser.add_argument('--seed-sigma',type=float,default=None,help='Opción de respaldo para sigma provisional')
    parser.add_argument('--seed',type=int,default=2026)
    parser.add_argument('--force',action='store_true',help='Reprocesar salidas IRAF existentes, guardando copia .bak')
    parser.add_argument('--timeout',type=int,default=3600,help='Tiempo máximo de cada ejecución IRAF')
    args=parser.parse_args();args.base=args.base.expanduser().resolve()
    if args.seed_fwhm<=0:parser.error('El FWHM inicial debe ser positivo.')
    if args.seed_sigma is not None and args.seed_sigma<=0:
        parser.error('El sigma inicial debe ser positivo.')
    if args.timeout<=0:parser.error('El timeout debe ser positivo.')
    if args.mode=='manual' and args.stage=='all':
        parser.error('Usa --stage bootstrap/select/thresholds/final/fixed/adaptive/compare individualmente en modo manual.')
    r=_cfgpaths(args.chip,args.filtro,args.base)
    stages=['bootstrap','select','thresholds','final','fixed','adaptive','compare'] if args.stage=='all' else [args.stage]
    for name in stages:
        print(f'\n====== {args.filtro} CHIP {args.chip} | {name.upper()} ======',flush=True)
        success=_stage(name,r,args)
        if success is False and args.mode=='manual':
            print('Se generó script IRAF. Ejecútalo y luego invoca la etapa siguiente.')
            break
    print('\nFinalizó la etapa solicitada.')


if __name__=='__main__':
    try: main()
    except (RuntimeError,ValueError,FileNotFoundError,subprocess.TimeoutExpired) as err:
        print(f'\nERROR: {err}',file=sys.stderr)
        sys.exit(2)
