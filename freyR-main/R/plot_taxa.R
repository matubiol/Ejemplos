#' Generates taxonomy bar plots and summarizes taxa relative abundances
#'
#' @importFrom phyloseq phyloseq
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq sample_sums
#' @importFrom phyloseq tax_table
#' @importFrom phyloseq psmelt
#' @importFrom phyloseq tax_glom
#' @importFrom phyloseq rank_names
#' @importFrom phyloseq transform_sample_counts
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr group_by_at
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_at
#' @importFrom dplyr mutate_if
#' @importFrom dplyr coalesce
#' @importFrom dplyr vars
#' @importFrom dplyr reframe
#' @importFrom dplyr filter
#' @importFrom plyr ddply
#' @importFrom tidyr spread
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 position_dodge2
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 %+%
#' @importFrom ggplot2 margin
#' @importFrom forcats fct_rev
#' @param data (Required). The list object imported with \code{\link{import_files}}.
#' @param var1 (Optional). Character string. Variable to be used by \code{\link{facet_wrap}}.
#'  Typically, the larger groups of samples.
#' @param var2 (Required). Character string. Samples will be grouped by this variable in the plot.
#'  Each bar in the plot will represent the elements of this variable.
#'  If you wish to represent individual samples, use \code{"Sample"}.
#'  Use the same as \code{"taxrank"} to make each taxon is represented by a single bar
#'  (e.g. \code{var2 = "Class", taxrank = "Class", type = "Taxa"}).
#' @param var3 (Required). Character string. Default \code{"Sample"}. A third variable only for taxa summary.
#'  Typically each individual sample.
#' @param type (Required). Character string. Default \code{"Sample"}. The plot type.
#'  The currently supported options are \code{c("Sample", "Taxa")}.
#'  Use the same as \code{"taxrank"} with \code{var2} to make each taxon is represented by a single bar
#'  (e.g. \code{var2 = "Class", taxrank = "Class", type = "Taxa"}).
#' @param taxrank (Required). Character string. Default \code{"Class"}. Taxonomic level to represent in the plot.
#' @param cutoff (Optional). Numeric.
#'  Determines which taxa will be grouped into ‘Others’ based on abundance (\code{top=FALSE}) or number of taxa (\code{top=TRUE}).
#' @param scales (Required). Character string. Default \code{"free_x"}. Customize axis scales.
#' @param top (Required). Logical. Default \code{FALSE}.
#'  If \code{TRUE}, the \code{cutoff} parameter determines the absolute number of displayed taxa + 'Others'.
#' @param norm (Required). Character string. Default \code{"rel"} .The normalization method.
#'  The currently supported options are \code{c("rel", "log10", "raw")}.
#' @param levels.1 (Optional). Character vector. Sort var1.
#' @param levels.2 (Optional). Character vector. Sort var2.
#' @param nrow (Required). Numeric. Default \code{1}.
#'  If \code{facet} is \code{TRUE}, the number of rows to be used by \code{\link{facet_wrap}}.
#' @param print (Required). Logical. Default \code{TRUE}. Print the plot.
#' @return Returns a list object and prints a plot.
#' @usage plot_taxa(data=NULL, var1=NULL, var2=NULL, var3="Sample", type=c("Sample","Taxa"), taxrank=c("Species","Genus",...),
#'  cutoff=NULL, scales="free_x", top=FALSE, norm=c("rel", "log10", "raw"), levels.1=NULL, levels.2=NULL, nrow=1, print=TRUE)
#' @examples
#' library(freyR)
#' # Import QIIME2 artifacts and metadata
#' files <- import_files(metadata="Path/metadata.tsv", table="Path/table.qza", taxonomy="Path/taxonomy.qza")
#'
#' # Plot taxonomy
#' taxa <- plot_taxa(files, var1="Location", var2="Treatments", levels.1=c("A", "B", "C"), levels.2=c("A", "B", "C"))
#' @export
plot_taxa <- function(data, var1=NULL, var2, var3="Sample", type="Sample", taxrank="Class", cutoff=NULL,
                      scales="free_x", top=FALSE, norm="rel", levels.1=NULL, levels.2=NULL, nrow=1, print=TRUE) {
  # Subset data components
  metadata <- data[["metadata"]]
  table <- data[["table"]]
  taxonomy <- data[["taxonomy"]]

  # Merge table, metadata and taxonomy into a phyloseq object
  rownames(metadata) <- NULL
  ps <- phyloseq(otu_table(table, taxa_are_rows = T),
                 tax_table(as.data.frame(taxonomy) %>% select(-ncol(taxonomy)) %>% column_to_rownames("Feature.ID") %>% as.matrix()),
                 sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))

  ps <-   prune_samples(sample_sums(ps) > 0, ps) #Remove empty samples
  taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (all(rank_names(ps) == taxa_levels)) {
    # How many features detected at each taxonomic level
    ps_df <- psmelt(ps) %>%
      filter(Abundance > 0) %>%
      distinct(OTU, .keep_all = TRUE) # De-replicate ASVs

    taxa_keys <- c("k__", "d__", "p__", "c__", "o__", "f__", "g__", "s__")

    raw <- sapply(ps_df[taxa_levels],
                  function(x) length(grep(paste0(taxa_keys, collapse = "|"), x)))

    perc <- sapply(ps_df[taxa_levels],
                   function(x) round(length(grep(paste0(taxa_keys, collapse = "|"), x))/length(taxonomy$Feature.ID)*100, 2))

    stats <- rbind(raw, perc)
    rownames(stats) <- c("Features", "%")
  } else {stats = NULL}

  # Normalization
  if (norm == "rel") {
    ylab = "Proportion"
    others_label <- paste0("Others (<", cutoff, "%)")
    cutoff1 <- cutoff/100
    ps1 <- transform_sample_counts(ps, function(x) x / sum(x))

  } else if (norm == "log10") {
    ylab = bquote(log[10]("No Reads"))
    others_label <- paste0("Others (<", cutoff, ")")
    cutoff1 <- cutoff
    ps1 <- transform_sample_counts(ps, function(x) log10(1+x))

  } else if (norm == "raw") {
    cutoff1 <- cutoff
    ylab = "No reads"
    others_label <- paste0("Others (<", cutoff, " reads)")
    ps1 = ps
  }

  # Agglomerate table to taxrank
  ps_df <- tax_glom(ps1, taxrank = taxrank, NArm = FALSE) %>%
    psmelt() %>%
    mutate_at(taxrank, function(x) {
      x <- gsub("Unknown", "Unassigned", x)
      x <- gsub(".+__", "", x)
      x <- ifelse(is.na(x), "Unassigned", x)
      return(x)
    })

  # Group data frame by taxrank, calculate average abundance
  means_table <- plyr::ddply(ps_df, as.formula(paste("~", taxrank)), function(x) c(mean=mean(x$Abundance)))
  means <- means_table
  ps1_df <- ps_df
  if (!is.null(cutoff)) {
    # Find taxa whose abundance is less than cutoff
    if (top == FALSE){
      others <- means[means$mean <= cutoff1,1]
    } else {
      others <- means[means$mean < arrange(means, desc(mean))[cutoff,2],1]
      others_label <- "Others"
    }

    # Change their name to "Others"
    means[is.element(means[,taxrank], others),taxrank] <- others_label
    ps1_df[is.element(as.data.frame(ps1_df)[,taxrank], others),taxrank] <- others_label
    ps1_df[,taxrank] <- factor(ps1_df[,taxrank], levels = c(means[order(-means$mean), taxrank][!(means[order(-means$mean), taxrank] %in% others_label)], others_label))
  } else {
    ps1_df[,taxrank] <- factor(ps1_df[,taxrank], levels = c(means[order(-means$mean), taxrank]))
  }

  ps1_df <- aggregate(as.formula(paste("Abundance ~", taxrank, "+ Sample")), ps1_df, FUN = sum) %>%
    merge(metadata, by.x = "Sample", by.y="SampleID")

  if (var2 == taxrank){
    ps2_df <- ps1_df %>%
      group_by_at(taxrank) %>%
      reframe(Abundance = mean(Abundance)) %>%
      as.data.frame()
  } else {
    ps2_df <- ps1_df %>%
      group_by_at(c(taxrank, var2)) %>%
      reframe(Abundance = mean(Abundance)) %>%
      as.data.frame()
  }

  # Make plot
  if (!is.null(levels.2)) {
    ps2_df[,var2] <- factor(ps2_df[,var2], levels = levels.2)
  }

  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#DDCC77", "#0072B2", "#D55E00", "#CC79A7",
                 "#332288", "#117733", "#44AA99", "#88CCEE", "#000000", "#CC6677", "#AA4499", "#882255")

  if (type == "Sample"){
    p <- ggplot(ps2_df, aes(x=.data[[var2]], y=Abundance, fill=.data[[taxrank]])) +
      geom_bar(stat = "identity") +
      labs(x = "", y = ylab) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(text=element_text(size=12),
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = margin(1,1,1,1.5, "cm"))
    
  } else if (type == "Taxa" & var2 == taxrank){
    p <- ggplot(ps2_df, aes(x=Abundance, y=fct_rev(.data[[taxrank]]), fill=.data[[taxrank]])) +
      geom_bar(stat = "identity", position=position_dodge2(width = 0.9, preserve = "single")) +
      labs(y = "", x = ylab) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, max(ps2_df$Abundance) * 1.1)) +
      theme(legend.position = "none")
    
  } else if (type == "Taxa" & var2 != taxrank){
    p <- ggplot(ps2_df, aes(x=Abundance, y=fct_rev(.data[[taxrank]]), fill=.data[[var2]])) +
      geom_bar(stat = "identity", position=position_dodge2(width = 0.9, preserve = "single")) +
      labs(y = "", x = ylab) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, max(ps2_df$Abundance) * 1.1))
  }
  p <- p  +
    scale_fill_manual(values=rep(cbPalette, times = 100)) +
    scale_color_manual(values=rep(cbPalette, times = 100))

  # Prepare tables
  means <- aggregate(as.formula(paste("mean ~", taxrank)), means_table, FUN = sum) %>% arrange(desc(mean))

  ps_df <- aggregate(as.formula(paste("Abundance ~", taxrank, "+ Sample")), ps_df, FUN = sum) %>%
    merge(metadata, by.x = "Sample", by.y="SampleID")

  means_v2 <- ps_df %>%
    group_by_at(c(taxrank, var2)) %>%
    reframe(Abundance = mean(Abundance)) %>%
    spread(var2, Abundance) %>%
    mutate_if(is.numeric, ~coalesce(., 0)) %>%
    as.data.frame()

  ps3_df <- ps_df %>%
    group_by_at(c(taxrank, var3)) %>%
    reframe(Abundance = mean(Abundance)) %>%
    as.data.frame()

  means_v3 <- ps3_df %>%
    select(c(Abundance, var3, taxrank)) %>%
    spread(var3, Abundance) %>%
    mutate_if(is.numeric, ~coalesce(., 0))

  taxa <- list(stats = stats, means = means, means_v2=means_v2, means_v3=means_v3)

  # Add var1
  if (!is.null(var1)) {
    ps2_df <- ps1_df %>%
      group_by_at(c(taxrank, var2, var1)) %>%
      reframe(Abundance = mean(Abundance)) %>%
      filter(Abundance > 0) %>%
      as.data.frame()

    means_v1 <- aggregate(as.formula(paste("Abundance ~", taxrank, "+", var1)), ps_df, FUN = mean) %>%
      spread(var1, Abundance) %>%
      mutate_if(is.numeric, ~coalesce(., 0))

    taxa[["means_v1"]] <- means_v1

    if (!is.null(levels.1)) {
      ps2_df[,var1] <- factor(ps2_df[,var1], levels = levels.1)
    }
    if (!is.null(levels.2)) {
      ps2_df[,var2] <- factor(ps2_df[,var2], levels = levels.2)
    }

    p <- p %+% ps2_df + facet_wrap(as.formula(paste("~",var1)), scales = scales, nrow = nrow)
    
    if (type == "Taxa") {
      p <- p + scale_x_continuous(expand = c(0, 0), limits = c(0, max(ps2_df$Abundance) * 1.1))
    }
  }
  taxa[["plot"]] <- p
  if (isTRUE(print)) {
    print(p)
  }
  return(taxa)
}
