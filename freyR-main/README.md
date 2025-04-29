# Instalación

```r
## Install Github packages
github_repo <- c("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", "jbisanz/qiime2R", "jfq3/QsRutils", "joey711/phyloseq")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

inst_github <- github_repo %in% installed.packages()
if(any(!inst_github)) {
    devtools::install_github(github_repo[!inst_github])
}

## Install freyR from gitlab
devtools::install_git("https://github.com/matubiol/Ejemplos/tree/main/freyR-main")
```

## Funciones principales del paquete `freyR`

### `import_files()`
Importa artefactos de QIIME 2 (`.qza`) y metadatos asociados, devolviendo una lista con las tablas de abundancia, taxonomía y metadatos.

**Ejemplo:**
```r
data <- import_files(metadata = "metadata.tsv", table = "table.qza", taxonomy = "taxonomy.qza")
```

### `filter_table()`
Filtra la tabla de características (`feature table`) según los metadatos proporcionados (por ejemplo, grupo de tratamiento), opcionalmente sobre la tabla normalizada.

**Ejemplo:**
```r
filtered <- filter_table(data, Treatment == "Control")
```

### `plot_alpha()`
Calcula y grafica métricas de diversidad alfa (e.g., Shannon, Simpson) en formato boxplot, incluyendo comparaciones estadísticas.

**Ejemplo:**
```r
plot_alpha(data, metric = "Shannon", x = "Treatment")
```

### `plot_beta()`
Calcula y visualiza la diversidad beta (e.g., PCoA, NMDS), incluyendo análisis como PERMANOVA o `pairwise.adonis`.

**Ejemplo:**
```r
plot_beta(data, method = "PCoA", dist = "bray", color = "Treatment")
```

### `plot_taxa()`
Genera gráficos de barras de composición taxonómica y resume las abundancias relativas a diferentes niveles taxonómicos.

**Ejemplo:**
```r
plot_taxa(data, rank = "Phylum", group = "Treatment")
```

### `plot_venn()`
Crea diagramas de Venn para comparar conjuntos de taxones entre grupos definidos por metadatos.

**Ejemplo:**
```r
plot_venn(data, group = "Treatment", rank = "Genus")
```
