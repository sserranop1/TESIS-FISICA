# Corrida nueva de Sextans A con aperturas por chip

Corrida completa del filtro V para los chips 2–8, con una carpeta de resultados separada. Cada subcarpeta `chips/chipN` contiene los productos de `sextans_auto.py`: selección de estrellas y cielo, controles, umbral DAOFIND, fotometría fija y adaptativa, comparación y resumen final.

| Chip | Radio de apertura |
|---:|---:|
| 2 | 2.5 × FWHM |
| 3 | 1.5 × FWHM |
| 4 | 1.5 × FWHM |
| 5 | 2 × FWHM |
| 6 | 6 × FWHM |
| 7 | 1 × FWHM |
| 8 | 1.5 × FWHM |

El anillo de cielo comienza en `max(5, múltiplo de apertura + 1) × FWHM`, para que no se superponga con la apertura. Por eso el chip 6 usa 7 × FWHM como radio interno del cielo; los demás usan 5 × FWHM.

`resumen_nueva_corrida.csv` reúne los parámetros principales y las métricas de comparación de esta corrida. Los CSV, catálogos, gráficas y resúmenes detallados están dentro de cada carpeta de chip. Los FITS V son enlaces a los originales en `dataSextansA/chips`, no copias.

El notebook de comparación por chip está en `../Comparacion_SextansA_corridas.ipynb`. Al ejecutarlo, regenera `comparacion_corridas_chip2_8.csv` y las figuras en `figuras_comparacion/`.

La instalación de IRAF tenía roto el enlace a su carpeta `extern`; para completar la corrida se usó `iraf_runtime/irafcl-wrapper`, con el entorno temporal de IRAF dentro de esta carpeta de resultados.

Para la comparación posterior con el procesamiento anterior, el resumen existente está en `../Resumen_SextansA_Filtro_V_Chips_1_a_8.csv`. Los resultados previos dentro de `dataSextansA` se conservaron.
