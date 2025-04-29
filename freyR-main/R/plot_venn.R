#' Generates Venn diagrams and provides list of taxa features
#'
#' @importFrom phyloseq phyloseq
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq tax_table
#' @importFrom phyloseq psmelt
#' @importFrom phyloseq tax_glom
#' @importFrom phyloseq prune_samples
#' @importFrom phyloseq prune_taxa
#' @importFrom phyloseq taxa_names
#' @importFrom phyloseq sample_names
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr filter
#' @importFrom ggVennDiagram ggVennDiagram
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 theme
#' @importFrom UpSetR upset
#' @importFrom UpSetR fromList
#' @importFrom writexl write_xlsx
#' @param data (Required). The list object imported with \code{\link{import_files}}.
#' @param norm (Required). Logical.
#' If \code{TRUE}, normalized table is used.
#' @param taxrank (Optional). Character string. Taxonomic level to represent in the diagram.
#' @param var1 (Required). Character string. Variable that determines the initial level of sample subsets.
#'  If \code{var2} is \code{NULL}, the number of circles in the Venn diagram will match the number of groups in \code{var1}.
#' @param var2 (Optional). Character string. Variable that determines the second level of sample subsets.
#' @param freq (Optional). Numeric.
#'  The percentage of the number of samples for each \code{var1/var2} in which a 'feature' must appear.
#' @param comb (Required). Logical. Default \code{FALSE}.
#'  Combine \code{var1} and \code{var2} in diagram. If the number of groups > 4, an \code{\link{upset}} plot will be created.
#' @param NArm (Optional). Logical. Default \code{TRUE}. CAUTION.
#'  The decision to prune (or not) taxa for which you lack categorical data could have a large effect on downstream analysis.
#' @param upset (Required). Logical. Default \code{FALSE}. Force an \code{\link{upset}} plot, even if the number of groups < 5.
#' @param names (Optional). Character string. This allows customization of the group names in the diagram.
#' @param sort (Optional). Character vector. Sort groups. If \code{names} is not \code{NULL}, use new names.
#' @param title (Optional). Character string. This add a title to the diagram.
#' @param write (Required). Logical. Default \code{FALSE}.
#' @param print (Required). Logical. Default \code{TRUE}. Print the plot.
#'  If \code{TRUE}, the list of features within each area of the diagram is written to an Excel file in the working directory.
#' @return Returns a list object and prints a diagram.
#' @usage plot_venn(data=NULL, norm=FALSE, taxrank=c("Species","Genus",...), var1=NULL, var2=NULL, freq=NULL,
#'  comb=FALSE, upset=FALSE, names=NULL, title=NULL, write=FALSE, print=TRUE)
#' @examples
#' library(freyR)
#' # Import QIIME2 artifacts and metadata
#' files <- import_files(metadata="Path/metadata.tsv", table="Path/table.qza", taxonomy="Path/taxonomy.qza")
#'
#' # Plot venn and get the features of each group.
#' venn <- plot_venn(files, taxrank = "Species", var1 = "Treatment", var2 = "Time", freq = 0.3, comb = FALSE)
#' @export
plot_venn <- function(data, norm, taxrank=NULL, var1, var2=NULL, freq=NULL, comb=FALSE, NArm=TRUE,
                      upset=FALSE, names=NULL, sort=NULL, title=NULL, write=FALSE, print=TRUE){
  #Subset data components
  metadata <- data[["metadata"]]
  if(isTRUE(norm)){
    table <- data[["table_norm"]]
  } else {
    table <- data[["table"]]
  }
  taxonomy <- data[["taxonomy"]]

  #Merge table, metadata and taxonomy into a phyloseq object
  rownames(metadata) <- NULL
  ps <- phyloseq(otu_table(table, taxa_are_rows = T),
                 tax_table(as.data.frame(taxonomy) %>% select(-ncol(taxonomy)) %>% column_to_rownames("Feature.ID") %>% as.matrix()),
                 sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))

  if (!is.null(taxrank)){
    ps1 <- tax_glom(ps, taxrank = taxrank, NArm = NArm)
  } else {ps1 = ps}

  metadata <- sample_data(ps1) %>% data.frame()

  # Filter features by frequency
  if (!is.null(freq)){
    feat <- c()
    if(is.null(var2)){
      for (x in unique(metadata[[var1]])){
        samples <- row.names(metadata)[metadata[var1] == x]
        psx <- prune_samples(sample_names(ps1) %in% samples, ps1)
        y <- phyloseq::filter_taxa(psx, function(z) sum(z > 1) > (freq*length(z)), TRUE) %>%
          otu_table() %>% row.names()
        feat <- c(feat, y)
      }
    } else {
      for (x in unique(metadata[[var1]])){
        samples.x <- row.names(metadata)[metadata[var1] == x]
        for (y in unique(metadata[[var2]])){
          samples.y <- row.names(metadata)[which(row.names(metadata)[metadata[var2] == y] %in% samples.x)]
          psx <- prune_samples(sample_names(ps1) %in% samples.y, ps1)
          feat.y <- phyloseq::filter_taxa(psx, function(z) sum(z > 1) > (freq*length(z)), TRUE) %>%
            otu_table() %>% row.names()
          feat <- c(feat, feat.y)
        }
      }
    }
    ps2 <- prune_taxa(taxa_names(ps1) %in% unique(feat), ps1)
  }

  ps_df <- ps1 %>%
    psmelt() %>%
    mutate(across(Kingdom:ncol(.), ~ifelse(is.na(.), "Unknown", .))) %>%
    filter(Abundance > 0)

  OTU_derep <- c()
  n=1
  for (i in ps_df[[taxrank]]){
    if (i != "Unknown"){
      OTU <- ps_df$OTU[ps_df[taxrank] == i][1]
      OTU_derep <- c(OTU_derep, OTU)
    } else {
      OTU <- ps_df$OTU[n]
      OTU_derep <- c(OTU_derep, OTU)
    }
    n=n+1
  }
  ps_df$OTU <- OTU_derep

  venn <- list()
  if(isFALSE(comb) | is.null(var2)){
    for (x in unique(ps_df[[var1]])){
      OTU <- ps_df$OTU[ps_df[,var1] == x]
      venn[[as.character(x)]] <- unique(OTU)
    }
  } else if (isTRUE(comb) & !is.null(var2)) {
    for (x in unique(ps_df[[var1]])){
      ps_df_x <- ps_df[ps_df[,var1] == x,]
      for(y in unique(ps_df_x[[var2]])){
        OTU <- ps_df_x$OTU[ps_df_x[,var2] == y]
        venn[[paste0(as.character(x), "-", as.character(y))]] <- unique(OTU)
      }
    }
  }
  if (!is.null(names)){
    names(venn) <- names
  }
  if (!is.null(sort)){
    venn <- venn[sort]
  }

  if (length(venn) > 4 | isTRUE(upset)) {
    p <- upset(fromList(venn), nsets = length(venn), order.by = "freq", nintersects = NA, sets = names(venn),
               keep.order = TRUE, main.bar.color = "#995ee1", sets.bar.color = "#995ee1", point.size = 2.5, line.size = 1,
               text.scale = c(1.5, 1.3, 1, 1, 1.5, 1.2))
  } else {
    p <- ggVennDiagram(venn) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
      ggtitle(title) + scale_fill_distiller(palette = "RdBu")
  }

  # Create a list of all the combinations
  combs <- unlist(lapply(1:length(venn), function(x) combn(names(venn), x, simplify = FALSE)), recursive = FALSE)

  names(combs) <- sapply(combs, function(x) paste0(x, collapse = "_"))

  Intersect <- function (x) {
    # Multiple set version of intersect
    # x is a list
    if (length(x) == 1) {
      unlist(x)
    } else if (length(x) == 2) {
      intersect(x[[1]], x[[2]])
    } else if (length(x) > 2){
      intersect(x[[1]], Intersect(x[-1]))
    }
  }

  Union <- function (x) {
    # Multiple set version of union
    # x is a list
    if (length(x) == 1) {
      unlist(x)
    } else if (length(x) == 2) {
      union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
      union(x[[1]], Union(x[-1]))
    }
  }

  Setdiff <- function (x, y) {
    # Remove the union of the y's from the common x's.
    # x and y are lists of characters.
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
  }

  elements <- lapply(combs, function(x) Setdiff(venn[x], venn[setdiff(names(venn), x)]))
  venn_data <- lapply(elements, function(x){
    y <- ps_df[ps_df$OTU %in% x, c(1, 3, which(colnames(ps_df) == "Kingdom"):ncol(ps_df))] %>%
      group_by(OTU) %>%
      mutate(Abundance = sum(Abundance)) %>%
      arrange(across(all_of(intersect(colnames(ps_df), rank_names(ps1))))) %>%
      distinct()
    names(y)[1] <- "FeatureID"
    unique(y)
  })
  if (isTRUE(write)){
    write_xlsx(venn_data, "venn_data.xlsx", format_headers = TRUE)
  }
  venn_data[["plot"]] <- p
  print(p)
  return(venn_data)
}
