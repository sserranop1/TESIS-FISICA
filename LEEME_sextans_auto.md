# Sextans A — Automatización de IRAF y fotometría de apertura

Script: `sextans_auto.py`. Basado en las funciones del notebook `procesamientoSextansA(20261026-155708).ipynb`.
No modifica el notebook original.

## Instalación (Mac, desde la Terminal normal, no desde `ecl>`)

Activa el entorno virtual con el que ejecutas el notebook de la tesis y comprueba que tiene `numpy`, `pandas`, `matplotlib`, `astropy` y `scipy`. Si falta alguno:

```bash
python -m pip install numpy pandas matplotlib astropy scipy
```

Comprueba que IRAF está instalado y accesible:

```bash
which irafcl
which ecl
```

Deja `sextans_auto.py` en `~/Desktop/Universidad/tesis/` o en cualquier carpeta. El script localiza por defecto `~/Desktop/Universidad/tesis/dataSextansA/chips/chipN/SexA_F_chipN.fits`.

## Ejecución automática completa de V · Chip 3

```bash
python sextans_auto.py --chip 3 --filtro V --mode auto --iraf-exe ecl
```

Si tu ejecutable es `irafcl`, usa `--iraf-exe irafcl`, o simplemente omite ese argumento.
En modo automático el programa emplea la opción `-f` del intérprete IRAF y registra un `.log` por etapa. **No se ha comprobado la invocación automática de IRAF en tu Mac**; si el intérprete no admite ese modo en tu instalación, utiliza el modo manual de abajo.

## Si prefieres seguir ejecutando los .cl en IRAF

```bash
python sextans_auto.py --chip 3 --filtro V --mode manual --stage bootstrap
```

Ejecuta el `.cl` que aparezca en la pantalla, dentro de `ecl>`; luego sigue:

```bash
python sextans_auto.py --chip 3 --filtro V --mode manual --stage select
python sextans_auto.py --chip 3 --filtro V --mode manual --stage thresholds
```

Ejecuta las siete pruebas de DAOFIND en `ecl>`, luego vuelve a ejecutar `--stage thresholds` **para leer los resultados** y calcular el threshold. Continúa:

```bash
python sextans_auto.py --chip 3 --filtro V --mode manual --stage final
python sextans_auto.py --chip 3 --filtro V --mode manual --stage fixed
python sextans_auto.py --chip 3 --filtro V --mode manual --stage adaptive
python sextans_auto.py --chip 3 --filtro V --mode manual --stage compare
```

Ejecuta el script `.cl` que genere cada etapa IRAF **antes** de continuar a la siguiente. Para `compare` no hace falta IRAF. En caso de que una ejecución IRAF esté pendiente, es normal que ese paso primero genere archivos `.cl`.

## Criterios de selección automática

- **Estrellas:** se inicia DAOFIND con FWHM preliminar de 4.5 píxeles (`--seed-fwhm` configurable), `threshold=5` y sigma provisional medido a partir de la imagen. Se requieren `0.2 <= SHARP <= 1.0`, `|SROUND| < 0.5`, `|GROUND| < 0.5`, separación a otras detecciones > 5 FWHM, sin relleno de 50000, pico <85% de SATURATE. Se selecciona una por cada celda de una malla de **5 columnas × 3 filas**, y se ajusta una gaussiana bidimensional. La **mediana de los 15 FWHM** será el FWHM operativo de IRAF.
- **Cielo:** se selecciona una caja de **15 × 15 píxeles** por celda (15 en total), lejos de detecciones y sin 50000. La estructura puntual local (*sharpness proxy* basado en residuo respecto a una imagen suavizada) debe ser baja; también se limitan los picos respecto al ruido, la proporción de píxeles altos y la dispersión. Se aplica el límite superior de Tukey a las 15 desviaciones estándar y se usa su mediana de aceptadas.
- **Fotometría:** conserva los criterios del notebook: `r_ap=3 FWHM`, `r_sky_in=5 FWHM`, ancho de anillo 15 píxeles, prueba de siete thresholds de 2 a 6, `MERR<=0.20`, aislamiento, `f_out` por Tukey y hasta cinco anillos separados 5 píxeles. Comprueba la reproducción exacta del control.

**Importante:** DAOFIND **no proporciona un SHARP para una región vacía**; su `sharpness` se define para objetos detectados. Por eso el script usa un indicador de estructura local bajo para el fondo, además de distancia a fuentes y pruebas de homogeneidad. El FWHM ajustado por una gaussiana **no es exactamente la columna `DIRECT` de `imexamine`**, así que revisa los valores y las figuras de control y, en una muestra, compáralos con tus mediciones manuales de los Chips 1–2. La selección de 15 objetos por malla es estricta: si una celda no tiene candidato aceptable, la etapa se detiene en lugar de inventar una medición.

## Resultados principales por chip

Se guardan en `chips/chipN` los archivos `<filtro>_chipN_auto_estrellas.csv`, `_auto_cielo.csv`, `_auto_control.png`, `_auto_config.json`, `_threshold.json`, catálogos de DAOFIND, archivos `_fixed.dat`, catálogos adaptativos, `_comparacion.csv`, `_movidas.csv`, `_contaminacion.csv` y `_resumen_final.json`.

Antes de confiar en un chip, **abre `_auto_control.png` y revisa los 15 círculos verdes y las 15 cajas azules**. Si las estrellas muestran saturación, perfiles extraños o las cajas cruzan nebulosidad/galaxias débiles, no aceptes esos parámetros automáticamente.

## Reprocesamiento

Si ya existen archivos, el modo automático los conserva y evita sobrescribirlos. Si necesitas rehacer una etapa, usa `--force` **solo después de verificar los datos**; el script conserva una copia `.bak` de los archivos de salida existentes.
