# Guion de presentación · 20 minutos

**Tesis:** Desarrollo de un software para fotometría de apertura con selección adaptativa del anillo de cielo  
**Autor:** Santiago Serrano Perico  
**Duración objetivo:** 20:00 · 30 diapositivas

El guion está escrito para hablar con calma mientras se señalan las figuras. Los tiempos incluyen pausas breves para dejar leer tablas y gráficos. Si hay que recuperar tiempo, resume la explicación de las diapositivas 4, 7 y 16; conserva el contexto metodológico de Sextans A y la interpretación cautelosa de la comparación.

## Distribución del tiempo

| Diapositivas | Bloque | Tiempo |
|---|---|---:|
| 1–5 | Apertura, pregunta, CCD y datos | 2:50 |
| 6–10 | Reducción, FWHM y detección | 3:25 |
| 11–20 | Fotometría, CMD y filtros de calidad | 6:20 |
| 21–25 | Estado de tesis y automatización Sextans A | 3:55 |
| 26–30 | Comparación, conclusiones y próximos pasos | 3:30 |
| **Total** |  | **20:00** |


## Diapositiva 1 · Desarrollo de un software para fotometría de apertura con selección adaptativa del anillo de cielo

**Tiempo:** 0:30

Buenos días. En esta presentación resumo el avance de mi tesis sobre fotometría de apertura. Primero mostraré el flujo que desarrollé y probé con M92 y M41: desde la reducción de imágenes hasta la construcción de una muestra fotométrica depurada. Después presentaré una extensión reciente aplicada a Sextans A, donde automatizo parte de la selección y comparo dos corridas con aperturas distintas por chip.

El objetivo de hoy es explicar tanto lo que ya está hecho como lo que falta validar.


## Diapositiva 2 · Una cadena: de la imagen cruda a una muestra fiable

**Tiempo:** 0:25

La historia del trabajo sigue esta cadena. Primero se corrigen las contribuciones del detector; luego se caracteriza el tamaño de las estrellas y el ruido del cielo. Con esos parámetros se detectan fuentes, se mide su flujo y se revisa la calidad de cada medición.

La tesis documenta sobre todo M92 y M41. Los archivos de procesamiento más recientes permiten añadir Sextans A como una extensión del flujo. Mantendré separadas esas dos etapas al interpretar los resultados.


## Diapositiva 3 · El reto: medir la estrella sin atribuirle el cielo

**Tiempo:** 0:50

En fotometría de apertura se suman las cuentas dentro de una región alrededor de la estrella y se resta una estimación del cielo local multiplicada por el área de la apertura. Por eso no basta con detectar una fuente: también hay que estimar correctamente el fondo.

Si otra estrella cae dentro del anillo de cielo, la estimación puede subir y el flujo neto cambia. Mi pregunta guía es cómo hacer reproducible esta medición y cómo identificar cuándo el anillo local deja de representar solamente el cielo. Esta pregunta motiva los dos filtros que mostraré más adelante.


## Diapositiva 4 · Del fotón al ADU: qué guarda cada píxel

**Tiempo:** 0:40

Un CCD convierte fotones incidentes en electrones acumulados. La electrónica lee esa carga y la convierte en un valor digital, normalmente expresado en ADU. Cada píxel contiene una combinación de señal estelar, luz del cielo, señal de lectura, corriente oscura y ruido.

La imagen cruda, entonces, no equivale directamente a una medida de brillo. Antes de hacer fotometría se corrigen los patrones instrumentales que pueden desplazar o modificar esas cuentas.


## Diapositiva 5 · Tres campos, dos etapas del trabajo

**Tiempo:** 0:25

Los primeros casos de estudio son M92, observado en B y V, y M41, cuyo archivo tenía varias extensiones que se separaron en canales blue, red y green. En ambos casos se trabajó con detección y fotometría de apertura.

La extensión reciente utiliza imágenes de Sextans A en el filtro V. La comparación nueva cubre los chips 2 a 8; el chip 1 no tiene una corrida nueva equivalente dentro de este conjunto.


## Diapositiva 6 · Antes de medir: corregir el detector

**Tiempo:** 0:40

Para M92 apliqué una secuencia de reducción CCD en IRAF. La imagen se revisa y se corrigen las contribuciones de lectura, se conserva la región útil, se resta el bias y se corrige la respuesta con el flat field.

El producto de esta etapa es una imagen más adecuada para localizar estrellas y medirlas. La reducción prepara los datos, pero no reemplaza la selección de parámetros fotométricos que depende de cada imagen.


## Diapositiva 7 · Cada corrección resuelve un efecto distinto del detector

**Tiempo:** 0:40

El overscan sirve para estimar el nivel electrónico de lectura. El trim retira los bordes que no pertenecen a la región útil. El bias corrige el patrón aditivo que aparece incluso en exposiciones de tiempo cero, y el flat field corrige diferencias relativas de sensibilidad entre píxeles.

Estas operaciones reducen efectos distintos, por eso se documentan como pasos separados. En las imágenes de M92 se aplicó esta preparación antes de pasar a la caracterización estelar y la fotometría.


## Diapositiva 8 · El FWHM conecta la forma estelar con el tamaño de apertura

**Tiempo:** 0:45

El FWHM es el ancho del perfil de una estrella a media altura y funciona como una escala del tamaño aparente de las fuentes en la imagen. En la tesis medí varias estrellas y utilicé la mediana para representar cada imagen.

La tabla muestra valores entre 2.22 y 3.36 píxeles para los ejemplos de M92 y M41. La apertura se fijó en tres veces esa mediana y se usó threshold 3.5 para la detección automática. Así los parámetros responden a la imagen y no a un radio fijo arbitrario.


## Diapositiva 9 · Elegir apertura es equilibrar flujo y ruido del cielo

**Tiempo:** 0:35

Una apertura demasiado pequeña puede dejar fuera parte de la luz de las alas estelares. Una apertura más grande recoge más flujo, pero también suma más píxeles de cielo y aumenta la posibilidad de incluir fuentes vecinas.

El criterio de tres veces el FWHM fue la elección base para M92 y M41. En Sextans A voy a mostrar una prueba más amplia: cambiar el múltiplo de apertura por chip para observar cómo cambia la muestra y la calidad fotométrica.


## Diapositiva 10 · El umbral se eligió estudiando cuántas fuentes detecta DAOFIND

**Tiempo:** 0:45

El umbral de DAOFIND controla cuánta señal debe sobresalir del ruido para registrar una fuente candidata. Al bajarlo, aparecen más detecciones, pero también aumenta la posibilidad de incluir falsos positivos. Al subirlo, se pierden fuentes débiles.

La curva resume las pruebas realizadas para estudiar ese compromiso. Para las imágenes trabajadas en la tesis se adoptó threshold 3.5, buscando recuperar fuentes sin que la muestra quede dominada por detecciones espurias.


## Diapositiva 11 · Las detecciones se inspeccionan antes de medirlas

**Tiempo:** 0:35

Aquí se ven detecciones superpuestas a imágenes de M92 y M41. DAOFIND entrega posiciones candidatas; tvmark permite mirar dónde quedaron las marcas y comprobar si corresponden a fuentes reales, bordes o regiones problemáticas.

Esta verificación visual es una parte importante del método: el catálogo automático se revisa antes de usarlo para calcular magnitudes. En M41 también fue necesario separar las extensiones y procesar cada canal por separado.


## Diapositiva 12 · De las coordenadas a la magnitud instrumental

**Tiempo:** 0:35

Una vez fijados los parámetros y revisadas las coordenadas, ejecuté phot o qphot para medir el flujo dentro de aperturas y estimar el fondo. Los resultados incluyen magnitud instrumental y error fotométrico para cada fuente.

En M92, las mediciones en B y V se emparejaron para construir el color B menos V. En M41, los canales separados permitieron probar el mismo tipo de flujo en un segundo campo.


## Diapositiva 13 · Las fuentes débiles tienen errores fotométricos mayores

**Tiempo:** 0:40

Estas dos gráficas muestran magnitud instrumental frente al error fotométrico para M92 y M41. En ambos casos, hacia magnitudes más débiles aumenta el error. Es el comportamiento esperado: las fuentes menos brillantes tienen menor señal respecto al ruido del cielo y del detector.

Por eso el catálogo no se interpreta mirando solo las magnitudes. También se revisa MERR, que ofrece una estimación formal de la incertidumbre de cada medición.


## Diapositiva 14 · MAG y MERR: brillo medido e incertidumbre estimada

**Tiempo:** 0:35

La magnitud instrumental se obtiene a partir del logaritmo del flujo neto medido. El error en magnitudes se relaciona con la razón señal a ruido; como aproximación, sigma magnitud es 1.0857 dividido por S sobre N.

MERR bajo representa menor incertidumbre formal bajo el modelo de medición. No demuestra por sí solo que la magnitud sea exacta frente a la realidad, porque para eso se necesita calibración o una referencia externa.


## Diapositiva 15 · M92: el primer diagrama color–magnitud

**Tiempo:** 0:45

Para construir el CMD instrumental se combinaron las mediciones de m92011 en V y m92015 en B. El eje horizontal es B menos V y el vertical es V. Cada punto representa una fuente con medidas emparejadas en ambos filtros.

Este primer diagrama permite explorar la distribución de color y brillo de la muestra, pero todavía no es un diagrama calibrado al sistema estándar. Es un resultado instrumental que depende de todo el procesamiento anterior.


## Diapositiva 16 · Cómo leer el CMD instrumental de M92

**Tiempo:** 0:35

En el eje horizontal, B menos V compara la magnitud de la fuente en dos bandas. En el vertical se representa V. Recordemos que la escala de magnitudes es inversa: un valor menor corresponde a una fuente más brillante.

El filtrado ayuda a distinguir mediciones potencialmente contaminadas y permite observar cómo se modifica la distribución. Por el carácter instrumental y parcial de la muestra, todavía no debemos inferir edades o parámetros físicos del cúmulo.


## Diapositiva 17 · Filtro 1: separar fuentes con vecinos demasiado cercanos

**Tiempo:** 0:40

El primer filtro calcula la distancia al vecino más cercano y la compara con tres veces el FWHM. Para m92011 el FWHM mediano es 3.355 píxeles, así que el límite es 10.065 píxeles. Se marcan como potencialmente afectadas las fuentes cuyo vecino está dentro de ese radio.

En esta muestra se rechazaron 14 estrellas por proximidad. La idea es reducir casos donde los perfiles estelares se mezclan y pueden sesgar una fotometría de apertura.


## Diapositiva 18 · Filtro 2: detectar fuentes dentro del anillo de cielo

**Tiempo:** 0:45

El segundo filtro mira el anillo utilizado para estimar el fondo. En el análisis de M92 se utilizó un radio interno de 15 píxeles y un ancho de 10, así que el anillo ocupa de 15 a 25 píxeles.

La figura muestra la fracción de píxeles altos por estrella. Las fuentes con valores atípicos en el anillo pueden tener una contribución de otra estrella o estructura. Con este criterio se marcaron 28 fuentes como contaminadas.


## Diapositiva 19 · Cómo se marca un anillo de cielo contaminado

**Tiempo:** 0:40

La regla de Tukey se aplica en dos niveles. Dentro de cada anillo, se considera alto un píxel que supera el tercer cuartil más 1.5 veces el rango intercuartílico. Luego se calcula la fracción de esos píxeles respecto al total del anillo.

Finalmente, la distribución de fracciones entre fuentes permite fijar un umbral para detectar anillos inusuales. Así el rechazo no depende de mirar cada caso a ojo, sino de una regla estadística reproducible que luego puede inspeccionarse visualmente.


## Diapositiva 20 · El CMD permite ver qué cambia después del filtrado

**Tiempo:** 0:30

El CMD coloreado muestra las mediciones que fueron marcadas por vecinos cercanos o por contaminación del anillo, junto con las restantes. Esto permite revisar si los filtros afectan ciertas regiones del diagrama de manera distinta.

Los dos grupos de rechazadas pueden superponerse; por eso los conteos de 14 y 28 no se deben sumar como si fueran objetos distintos. El resultado es una muestra mejor caracterizada, no una garantía de calibración.


## Diapositiva 21 · Estado de la tesis: la base fotométrica está construida

**Tiempo:** 0:40

En el manuscrito ya están documentadas la reducción de M92, la separación de canales de M41, la caracterización FWHM y cielo, la detección, la fotometría, el CMD y los dos criterios de calidad.

El trabajo sigue en desarrollo. Entre las etapas pendientes están calibrar magnitudes, interpretar un CMD calibrado, comparar con fotometría PSF en regiones densas e integrar al manuscrito los resultados recientes de Sextans A.


## Diapositiva 22 · Extensión reciente: selección automática para Sextans A

**Tiempo:** 0:50

Los archivos recientes aplican el flujo a imágenes V de Sextans A por chip. El script organiza la selección de parámetros, la detección, la fotometría fija y el análisis de contaminación del cielo. También guarda tablas y figuras de control para revisar cada etapa.

En esta imagen del chip 3 se marcan quince estrellas usadas para estimar el FWHM y quince cajas para caracterizar el cielo. La revisión visual ayuda a detectar si la automatización eligió regiones razonables antes de confiar en los parámetros.


## Diapositiva 23 · La selección automática se apoya en controles distribuidos por el chip

**Tiempo:** 0:50

Para el FWHM se parte de candidatas estelares con forma adecuada, distribuidas por una grilla de cinco columnas y tres filas. A cada una se le ajusta un perfil gaussiano bidimensional y se rechazan ajustes atípicos antes de tomar la mediana.

Para el cielo se evalúan quince cajas de 15 por 15 píxeles alejadas de fuentes. También se filtran valores atípicos y se usa la mediana aceptada de sigma. El gráfico de control registra tanto las posiciones como los parámetros derivados.


## Diapositiva 24 · La nueva corrida usa un múltiplo de FWHM distinto por chip

**Tiempo:** 0:45

La solicitud fue probar aperturas de 2.5 FWHM para chip 2, 1.5 para chips 3 y 4, 2 para chip 5, 6 para chip 6, 1 para chip 7 y 1.5 para chip 8. En el gráfico también aparecen los radios en píxeles derivados del FWHM de cada chip.

El chip 6 alcanza 27.72 píxeles de radio de apertura. El chip 1 queda fuera porque no se procesó en esta nueva corrida comparable.


## Diapositiva 25 · El anillo adaptativo prueba posiciones sucesivas del cielo

**Tiempo:** 0:50

El anillo parte después de la apertura. En general, el radio interno inicia en cinco FWHM, pero si la apertura es grande se desplaza hasta al menos un FWHM más allá del radio de apertura. El ancho del anillo es de quince píxeles.

Para cada estrella se evalúan hasta cinco posiciones, moviéndose cinco píxeles hacia afuera. Se escoge la primera que pasa el límite de contaminación. En chip 6, con apertura de seis FWHM, el cielo empieza en siete FWHM para que las regiones no se crucen.


## Diapositiva 26 · Qué cambia entre la corrida anterior y la nueva

**Tiempo:** 0:40

La comparación usa solo chips 2 a 8 en V. La corrida anterior empleó tres FWHM de apertura en todos estos chips. Para chips 3 a 5 y 7 a 8, la nueva prueba modifica el múltiplo de apertura.

Hay dos salvedades: chip 2 cambia de una corrida manual en IRAF, con seis estrellas y 22 cajas, a una automática con 15 estrellas y 15 cajas. Chip 6 además desplaza el anillo interno de cinco a siete FWHM. Así que comparamos resultados de corridas completas, no un único cambio aislado.


## Diapositiva 27 · La cantidad de mediciones válidas cambió de forma distinta por chip

**Tiempo:** 0:45

El gráfico muestra el cambio porcentual y los conteos anterior y nuevo de fotometrías válidas. En total, para chips 2 a 8, se pasa de 19,558 a 18,880 registros, una disminución cercana al 3.5 por ciento.

El chip 2 aumenta 6 por ciento; los otros seis chips disminuyen entre aproximadamente 1 y 8 por ciento. Estos son conteos de muestra, no una medida directa de precisión. La diferencia de método en chip 2 también debe tenerse en cuenta.


## Diapositiva 28 · La calidad no mejora igual en todos los chips

**Tiempo:** 0:50

Aquí se comparan dos métricas entre corridas para fuentes con anillo desplazado y fotometría válida. La mediana de la diferencia absoluta entre magnitud adaptativa y fija disminuye o empata en los siete chips. Sin embargo, el MERR adaptativo baja solamente en los chips 3, 4 y 7; sube en 2, 5, 6 y 8.

No hay una respuesta universal por chip. La selección de apertura cambia la muestra y el fondo local, así que los indicadores deben leerse juntos y para su subconjunto correspondiente.


## Diapositiva 29 · Qué sabemos — y qué todavía no concluye la comparación

**Tiempo:** 0:35

La comparación muestra que el radio de apertura cambia el número de mediciones y que el desplazamiento adaptativo del cielo puede cambiar magnitudes para fuentes cuyo anillo inicial está contaminado. También muestra que las medianas de error no mejoran de manera uniforme.

Sin magnitudes estándar o un catálogo de referencia, no podemos afirmar qué corrida es más exacta. El siguiente paso es validar el comportamiento con datos externos o experimentos controlados y registrar esos límites en la tesis.


## Diapositiva 30 · Una base reproducible; ahora falta calibrarla y validarla

**Tiempo:** 0:40

Como conclusión, el trabajo ya conecta imágenes CCD con una muestra fotométrica depurada. M92 y M41 permitieron construir y documentar el flujo base. Sextans A añade automatización por chip y una comparación sistemática de radios de apertura y anillos de cielo.

Los próximos pasos son integrar este análisis al manuscrito, calibrar magnitudes con punto cero y correcciones atmosféricas y validar contra fuentes de referencia. En campos densos también quiero comparar apertura con PSF. Gracias.


## Referencias visuales y fuentes

Las figuras de M92, M41 y sus filtros provienen de los archivos de la tesis. Los gráficos comparativos de Sextans A se construyeron a partir de los CSV de ambas corridas y del cuaderno `Comparacion_SextansA_corridas.ipynb`. La imagen de Sextans A se acredita a Subaru Telescope / NAOJ.
