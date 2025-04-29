# Ejemplos de scripts implementando diferentes lenguajes de programación

Este repositorio contiene ejemplos de scripts desarrollados en diferentes lenguajes de programación para el análisis de datos biológicos.

## 1. `denovo_assembly_ONT.sh` (Bash)

Este script realiza un ensamblaje **de novo** a partir de datos de secuenciación de **Oxford Nanopore**, siguiendo varias etapas de procesamiento:

- **Preprocesamiento de lecturas:** Copia las lecturas a una carpeta de trabajo, elimina adaptadores con **Porechop**, filtra por calidad con **Chopper**, y genera estadísticas de calidad con **NanoPlot** y **SeqKit**.
- **Ensamblaje:** Usa **Flye** para ensamblar las lecturas y **Medaka** para pulir el ensamblaje.
- **Mapeo y evaluación:** Mapea las lecturas contra el ensamblaje con **Minimap2** y genera estadísticas de alineamiento con **Samtools**. Finalmente, usa **QUAST** para evaluar la calidad del ensamblaje.
- **Pulido opcional con Illumina:** Si se proporciona un directorio con lecturas de **Illumina**, el script realiza una fase adicional de pulido con **Pilon** para mejorar la precisión del ensamblaje.

## 2. `taxa_check.py` (Python)

Este script evalúa la **disponibilidad y resolución taxonómica** de secuencias en bases de datos de referencia. Sus principales funciones incluyen:

- **Carga y procesamiento de taxonomías de referencia**, incluyendo la extracción de información desde bases de datos como **NCBI** y bases de datos internas.
- **Búsqueda de secuencias** correspondientes a un taxón específico en la base de datos de referencia, con la opción de descargar secuencias de **NCBI** si no están presentes.
- **Generación de reportes** sobre la cantidad de secuencias encontradas, su nivel taxonómico más detallado y su origen (**NCBI** o una base de datos específica).
- **Clasificación taxonómica de secuencias** usando diferentes métodos de asignación (**Naïve Bayes** o **VSEARCH**) en **QIIME 2**.
- **Creación de archivos de salida** en formatos tabulados y **Excel** con información sobre la precisión taxonómica alcanzada para cada secuencia analizada.

## 3. `freyR-main` (R)

Es un **paquete en R** diseñado para el análisis y visualización de datos de **microbiomas**, proporcionando herramientas para la **manipulación, exploración y representación gráfica** de datos de diversidad microbiana a partir de conteos taxonómicos. Facilita la integración de datos **filogenéticos, taxonómicos y ambientales**, permitiendo **análisis estadísticos y comparaciones entre muestras**, similar a **phyloseq**.

