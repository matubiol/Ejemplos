# Ejemplos de scripts implementando diferentes lenguajes de programación

Este repositorio contiene ejemplos de scripts desarrollados en diferentes lenguajes de programación para el análisis de datos biológicos.

## 1. `denovo_assembly_ONT.sh` (Bash)

Este script realiza un ensamblaje **de novo** a partir de datos de secuenciación de **Oxford Nanopore**, siguiendo varias etapas de procesamiento:

- **Preprocesamiento de lecturas:** Copia las lecturas a una carpeta de trabajo, elimina adaptadores con **Porechop**, filtra por calidad con **Chopper**, y genera estadísticas de calidad con **NanoPlot** y **SeqKit**.
- **Ensamblaje:** Usa **Flye** para ensamblar las lecturas y **Medaka** para pulir el ensamblaje.
- **Mapeo y evaluación:** Mapea las lecturas contra el ensamblaje con **Minimap2** y genera estadísticas de alineamiento con **Samtools**. Finalmente, usa **QUAST** para evaluar la calidad del ensamblaje.
- **Pulido opcional con Illumina:** Si se proporciona un directorio con lecturas de **Illumina**, el script realiza una fase adicional de pulido con **Pilon** para mejorar la precisión del ensamblaje.

### Ejecución básica

```bash
bash denovo_assembly_ONT.sh --dir_reads=/ruta/a/lecturas_ONT --output=/ruta/salida
```
### Ejecución con pulido adicional usando lecturas Illumina
```bash
bash denovo_assembly_ONT.sh --dir_reads=/ruta/a/lecturas_ONT --output=/ruta/salida --illumina=/ruta/a/lecturas_illumina
```
**Parámetros**

--dir_reads=: Ruta a la carpeta que contiene archivos .fastq.gz de lecturas Nanopore.

--output=: Ruta de la carpeta de salida donde se guardarán los resultados.

--illumina=: (Opcional) Ruta a la carpeta con lecturas Illumina emparejadas (*_R1_*.fastq.gz, *_R2_*.fastq.gz) para el pulido con Pilon.

## 2. `taxa_check.py` (Python)

Este script evalúa la **disponibilidad y resolución taxonómica** de secuencias en bases de datos de referencia. Sus principales funciones incluyen:

- **Carga y procesamiento de taxonomías de referencia**, incluyendo la extracción de información desde bases de datos como **NCBI** y bases de datos internas.
- **Búsqueda de secuencias** correspondientes a un taxón específico en la base de datos de referencia, con la opción de descargar secuencias de **NCBI** si no están presentes.
- **Generación de reportes** sobre la cantidad de secuencias encontradas, su nivel taxonómico más detallado y su origen (**NCBI** o una base de datos específica).
- **Clasificación taxonómica de secuencias** usando diferentes métodos de asignación (**Naïve Bayes** o **VSEARCH**) en **QIIME 2**.
- **Creación de archivos de salida** en formatos tabulados y **Excel** con información sobre la precisión taxonómica alcanzada para cada secuencia analizada.

### Ejecución mínima (con parámetros explícitos)

```bash
python taxa_check.py \
  --taxon "Escherichia coli" \
  --amplicon 16S \
  --trimmed_db Y \
  --classifier path/a/classifier.qza \
  --ref_seqs_path path/a/secuencias.qza \
  --ref_taxa_path path/a/taxonomia.qza \
  --primer_f ACTCCTACGGGAGGCAGCAG \
  --primer_r GGACTACHVGGGTWTCTAAT \
  --thread 8 \
  --taxa_approach nb
```
### Ejecución usando archivo de configuración

```bash
python taxa_check.py \
  --taxon "Homo sapiens" \
  --amplicon 16S \
  --trimmed_db N \
  --config config_file.tsv \
  --thread 8 \
  --taxa_approach both \
  --download_all Y
```

**Parámetros**

--taxon: Nombre del taxón de interés. Puede ser un nombre o un archivo .txt con una lista de taxones (uno por línea). Obligatorio.

--amplicon: Región amplicónica objetivo (por ejemplo, 16S, ITS, COI, 18S). Obligatorio.

--trimmed_db: Indica si la base de datos ya está recortada (Y o N). Obligatorio.

--classifier: Ruta al clasificador .qza (necesario si no se usa --config).

--ref_seqs_path: Ruta a las secuencias de referencia (.qza) (necesario si no se usa --config).

--ref_taxa_path: Ruta a la taxonomía de referencia (.qza) (necesario si no se usa --config).

--primer_f: Secuencia del cebador forward (necesario si no se usa --config).

--primer_r: Secuencia del cebador reverse (necesario si no se usa --config).

--config: Archivo .tsv con información predefinida para cada amplicón (en lugar de especificar todos los parámetros anteriores).

--thread: Número de hilos/CPU a usar (por defecto: 50% de los disponibles).

--taxa_approach: Método de clasificación taxonómica: nb (Naïve Bayes), vs (VSEARCH) o both. Por defecto: both.

--download_all: Si no se encuentra el taxón en la base de datos de referencia, indica si se deben descargar automáticamente secuencias desde NCBI (Y o N). Por defecto: N.

## 3. `freyR-main` (R)

Es un **paquete en R** diseñado para el análisis y visualización de datos de **microbiomas**, proporcionando herramientas para la **manipulación, exploración y representación gráfica** de datos de diversidad microbiana a partir de conteos taxonómicos. Facilita la integración de datos **filogenéticos y taxonómicos**, permitiendo **análisis estadísticos y comparaciones entre muestras**, similar a **phyloseq**.

