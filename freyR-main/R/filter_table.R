#' Filter feature table based on metadata
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @param data (Required). The list object imported with \code{\link{import_files}}.
#' @param norm (Required). Logical. Default \code{FALSE}.
#' If \code{TRUE}, filter normalized table.
#' @param ... (Required). The sub-setting expression that can be passed to the function.
#'  These arguments may be used to filter samples based on specific metadata.
#'  For example, \code{Treatment == "Group A"} can be used to filter samples where the treatment is equal to "Group A".
#' @return List with the filtered table.
#' @usage subsetSamples(data=NULL, ..., norm=FALSE)
#' @examples
#' library(freyR)
#' # Import QIIME2 artifacts and metadata
#' files <- import_files(table="Path/table.qza", metadata="Path/metadata.tsv")
#'
#' # Filter table based on metadata
#' table_filtered <- subsetSamples(files, Treatment == "Group A", norm = FALSE)
#' @export
subsetSamples <- function(data, ..., norm){
  #Subset data components
  metadata <- data[["metadata"]]
  if(isTRUE(norm)){
    table <- data[["table_norm"]]
    n = which(names(data) == "table_norm")
  } else {
    table <- data[["table"]]
    n = which(names(data) == "table")
  }
  sample_ids <- as.character(subset(metadata, ...)$SampleID)
  table_filt <- table %>% select(any_of(sample_ids))
  data[[n]] <- table_filt

  return(data)
}

#' Filter feature table based on taxonomy
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate_all
#' @param data (Required). The list object imported with \code{\link{import_files}}.
#' @param ... (Optional). If \code{by = "Taxa"}, the sub-setting expression that can be passed to the function.
#'  These arguments may be used to filter samples based on taxonomy.
#'  For example, Genus == "Escherichia" can be used to filter features where the genus is equal to "Escherichia".
#' @param freq (Optional). If \code{by = "Freq"}, a function or a list of functions that take a vector of
#'  abundance values and return a logical value. For example, function(x) sum(x > 1) > (0.3*length(x)) can be used
#'  to filter features seen in at least 30 \% of the samples.
#' @param by (Required). Character string. Filter features based on taxonomy or frequency.
#'  The currently supported options are \code{c("Taxa", "Freq")}.
#' @param norm (Required). Logical.
#'  If \code{TRUE}, filter normalized table.
#' @return List with the filtered table.
#' @usage filter_samples(data=NULL, ..., freq = NULL, by=c("Taxa","Freq"), norm=FALSE)
#' @examples
#' library(freyR)
#' # Import QIIME2 artifacts and metadata
#' files <- import_files(table="Path/table.qza", taxonomy="Path/taxonomy.qza")
#'
#' # Filter table based on taxonomy
#' table_filtered <- subsetTaxa(files, Genus == "Escherichia", by="Taxa", norm = FALSE)
#'
#' # WARNING!
#' # When filtering out non-target taxa, this function will also remove any NA values by default.
#' # To include NA values in the subset, you need to modify your filter condition to explicitly allow NA values.
#' # This ensures that taxa with NA in the specified level are retained during the filtering process:
#' table_filtered <- filter_taxa(files_target, !Species == non_target | is.na(Species), by = "Taxa", norm = FALSE)
#'
#' # Filter table based on frequency
#' #Remove taxa not seen more than 3 times in at least 20% of the samples.
#' table_filtered <- subsetTaxa(files, freq = function(x) sum(x > 1) > (0.3*length(x)), by="Freq", norm = FALSE)
#' @export
subsetTaxa <- function(data, ..., freq, by, norm){
  #Subset data components
  taxonomy <- data[["taxonomy"]]
  if(isTRUE(norm)){
    table <- data[["table_norm"]]
    n = which(names(data) == "table_norm")
  } else {
    table <- data[["table"]]
    n = which(names(data) == "table")
  }
  if (by == "Taxa"){
    taxonomy <- taxonomy %>%
      mutate_all(as.factor) %>%
      mutate_all(~gsub(".+__", "", .))
    taxa_ids <- as.character(subset(taxonomy, ...)$Feature.ID)
    table_filt <- table %>% dplyr::filter(rownames(.) %in% taxa_ids)
    data[[n]] <- table_filt
  } else if (by == "Freq"){
    feat <- apply(table, 1, freq)
    table_filt <- table[feat, ]
    data[[n]] <- table_filt
  }
  return(data)
}
