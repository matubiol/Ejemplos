#' Import QIIME2 artifacts and metadata
#'
#' @importFrom qiime2R read_q2metadata
#' @importFrom qiime2R read_qza
#' @importFrom tidyr as_tibble
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @param metadata (Required). A path to the Q2-formatted metadata.
#' @param table (Optional). A path to the Q2-feature-table artifact.
#' @param table_norm (Optional). A path to the normalized Q2-feature-table artifact.
#' @param taxonomy (Optional). A path to the Q2-taxonomy artifact for features in the provided feature table(s).
#' @return A list object with the imported files.
#' @usage import_files(metadata="Path/metadata.tsv", table=NULL, table_norm=NULL, taxonomy=NULL)
#' @examples
#' library(freyR)
#' # Import QIIME2 artifacts and metadata
#' files <- import_files(metadata="Path/metadata.tsv",
#'                       table="Path/table_norm.qza",
#'                       table_norm="Path/table_norm.qza",
#'                       taxonomy="Path/taxonomy.qza")
#' @export
import_files <- function(metadata, table=NULL, table_norm=NULL, taxonomy=NULL){
  cat("Starting import_files...\n")

  files <- list(metadata = read_q2metadata(metadata))
  cat("Metadata loaded successfully.\n")

  if (!is.null(table)){
    table <- read_qza(table)
    table = as.data.frame(table$data)
    files[["table"]] <- table
    cat("Table loaded successfully.\n")
  }

  if (!is.null(table_norm)){
    table_norm <- read_qza(table_norm)
    table_norm = as.data.frame(table_norm$data)
    files[["table_norm"]] <- table_norm
    cat("Normalized table loaded successfully.\n")
  }

  if (!is.null(taxonomy)){
    taxonomy <- read_qza(taxonomy)
    taxonomy <- taxonomy$data %>%
      as_tibble() %>%
      separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
    files[["taxonomy"]] <- taxonomy
    cat("Taxonomy loaded successfully.\n")
  }
  return(files)
}
