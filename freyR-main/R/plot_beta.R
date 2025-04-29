#' Calculate and plot beta diversity results
#'
#' @importFrom phyloseq phyloseq
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq plot_ordination
#' @importFrom magrittr %>%
#' @importFrom tibble column_to_rownames
#' @importFrom vegan vegdist
#' @importFrom vegan adonis2
#' @importFrom pairwiseAdonis pairwise.adonis
#' @importFrom ape pcoa
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 geom_abline
#' @importFrom plotly ggplotly
#' @importFrom cowplot ggdraw
#' @importFrom cowplot get_legend
#' @importFrom cowplot plot_grid
#' @importFrom cowplot draw_label
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid nullGrob
#' @param data (Required). The list object imported with \code{\link{import_files}}.
#' @param var1 (Required). Character string. Main variable to evaluate.
#'  Stats results in the plot reference to this variable. This variable determines \code{color} in the plot.
#'  If the number of groups within this variable exceeds the length of the palette,
#'  a combination of colors and shapes will be used to represent the variable.
#'  Consequently, it will not be possible to represent a second variable (\code{var2}) using shapes.
#' @param var2 (Optional). Character string. This variable determines \code{shape} in the plot.
#' @param weighted (Required). Logical. Default \code{FALSE}. Use weighted distance. This generates a second plot.
#' @param stats (Required). Logical. Default \code{TRUE}. Include PERMANOVA results in the plot.
#' @param dist.u (Required). Character string. Default \code{"jaccard"}.
#'  The dissimilarity index to be used in the binary version. See \code{\link{vegdist}} for more information.
#' @param dist.w (Required). Character string. Default \code{"jaccard"}.
#'  The dissimilarity index to be used in the non-binary version. See \code{\link{vegdist}} for more information.
#' @param label (Optional). Character string.
#'  Use a third variable represented by \code{label} in the plot.
#' @param inter (Required). Logical. Default \code{FALSE}.
#'  If \code{var2} is not \code{NULL}, \code{inter} determines if the \code{\link{adonis2}} formula
#'  will be set with interaction of factors (*) or not (+).
#' @param title (Optional). Character vector. Title for each plot.
#' @param print (Required). Logical. Default \code{TRUE}. Print the plot.
#' @return Returns a list object and prints a plot.
#' @usage plot_beta(data=NULL, var1=NULL, var2=NULL, weighted=FALSE, stats=TRUE,
#'  dist.u="jaccard", dist.w="jaccard", label=NULL, inter=FALSE, print=TRUE)
#' @examples
#' library(freyR)
#' # Import QIIME2 artifacts and metadata
#' files <- import_files(metadata="Path/metadata.tsv", table_norm="Path/table_norm.qza", taxonomy="Path/taxonomy.qza")
#'
#' # Calculate and plot beta diversity results
#' beta <- plot_beta(files, var1="Treatment", var2="Location", title = c("uJaccard - Bacteria", "wJaccard - Bacteria"))
#' @export
plot_beta <- function(data, var1, var2=NULL, dist.u="jaccard", weighted=FALSE, dist.w="jaccard",
                      label=NULL, stats=TRUE, inter=FALSE, title=c(" ", " "), print=TRUE) {
  #Subset data components
  metadata <- data[["metadata"]]
  table_norm <- data[["table_norm"]]
  
  #Merge table and metadata into a phyloseq object
  table_norm <- table_norm[colnames(table_norm) %in% metadata$SampleID[metadata[var1] != ""]]
  rownames(metadata) <- NULL
  ps <- phyloseq(otu_table(table_norm, taxa_are_rows = T),
                 sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))
  
  if (var1 == "Sample") {
    stats = FALSE
    ps@sam_data[["Sample"]] <- rownames(data.frame(ps@sam_data))
  }
  
  #Unweighted distances and PCoA ordination
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                 "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")
  
  if (length(unique(ps@sam_data[["Location"]])) > length(cbPalette)) {
    shape = var1
    scale_shape <- scale_shape_manual(values = rep(c(15:18), each = length(cbPalette)))
    guides <- guides(color = guide_legend(override.aes = aes(label = "")))
  } else {
    shape = var2
    scale_shape <- NULL
    guides <- guides(color = guide_legend(order = 1, override.aes = aes(label = "")),
                     shape = guide_legend(order = 2))
  }
  
  dis_uw <- vegdist(t(otu_table(ps)), dist.u, binary = TRUE)
  PCoA_unweighted <- pcoa(dis_uw)
  
  p.u <- plot_ordination(ps, PCoA_unweighted, color = var1, shape = shape, label = label) +
    scale_colour_manual(values=rep(cbPalette, times = 4)) +
    scale_shape +
    guides +
    geom_point(size=2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_light()
  
  if (isTRUE(weighted)) {
    dis_w <- vegdist(t(otu_table(ps)), dist.w, binary = FALSE)
    PCoA_weighted <- pcoa(dis_w)
    p.w <- plot_ordination(ps, PCoA_weighted, color = var1, shape = shape, label = label) +
      scale_colour_manual(values=rep(cbPalette, times = 4)) +
      scale_shape +
      guides +
      geom_point(size=2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      theme_light()
    
    legend <- suppressWarnings(get_legend(p.u))
    grouped_plots <- plot_grid(p.u + labs(title = title[1]) + theme(legend.position="none"),
                               p.w + labs(title = title[2]) + theme(legend.position="none"),
                               ncol = 2, align = "hv")
    p <- plot_grid(grouped_plots, legend, rel_widths = c(.4, .05))
  } else {
    p <- p.u + labs(title = title[1])
  }
  if (stats == TRUE) {
    if(!is.null(var2) && inter == TRUE){
      formula <- paste(var1, "*", var2)
      row <- 3
    } else if (!is.null(var2)){
      formula <- paste(var1, "+", var2)
      row <- 1
    } else {
      formula <- paste(var1)
      stats.v2 <- NULL
      row <- 1
    }
    u.adonis <- adonis2(as.formula(paste("dis_uw ~", formula)), data = data.frame(sample_data(ps)))
    u.pairwiseAdonis_v1 <- pairwise.adonis(dis_uw, ps@sam_data[[var1]], p.adjust.m = "BH")
    beta_stats <- list(u.adonis = as.data.frame(u.adonis),
                       u.pairwiseAdonis_v1 = cbind(data.frame(Pairs = u.pairwiseAdonis_v1$pairs),
                                                   round(u.pairwiseAdonis_v1[,5:7], 3)))
    if (!is.null(var2)){
      u.pairwiseAdonis_v2 <- pairwise.adonis(dis_uw, ps@sam_data[[var2]], p.adjust.m = "BH")
      beta_stats[["u.pairwiseAdonis_v2"]] <- cbind(data.frame(Pairs = u.pairwiseAdonis_v2$pairs),
                                                   round(u.pairwiseAdonis_v2[,5:7], 3))
    }
    if (weighted == TRUE) {
      w.adonis <- adonis2(as.formula(paste("dis_w ~", formula)), data = data.frame(sample_data(ps)))
      w.pairwiseAdonis_v1 <- pairwise.adonis(dis_w, ps@sam_data[[var1]], p.adjust.m = "BH")
      beta_stats[["w.adonis"]] <- as.data.frame(w.adonis)
      beta_stats[["w.pairwiseAdonis_v1"]] <- cbind(data.frame(Pairs = w.pairwiseAdonis_v1$pairs),
                                                   round(w.pairwiseAdonis_v1[,5:7], 3))
      if (!is.null(var2)){
        w.pairwiseAdonis_v2 <- pairwise.adonis(dis_uw, ps@sam_data[[var2]], p.adjust.m = "BH")
        beta_stats[["w.pairwiseAdonis_v2"]] <- cbind(data.frame(Pairs = w.pairwiseAdonis_v2$pairs),
                                                     round(w.pairwiseAdonis_v2[,5:7], 3))
      }
      p.u2 <- ggdraw(p.u + theme(legend.position="none") + labs(title = title[1])) +
        draw_label(bquote(paste(p == .(round(u.adonis[row,"Pr(>F)"],3)),"; ", R^2 == .(round(u.adonis[row,"R2"],3)))),
                   x = Inf, y = Inf, hjust = 1.05, vjust = 1.15, size = 14)
      
      p.w2 <- ggdraw(p.w + theme(legend.position="none") + labs(title = title[1])) +
        draw_label(bquote(paste(p == .(round(w.adonis[row,"Pr(>F)"],3)),"; ", R^2 == .(round(w.adonis[row,"R2"],3)))),
                   x = Inf, y = Inf, hjust = 1.05, vjust = 1.15, size = 14)
      
      if(!is.null(var2) && inter == FALSE){
        p.u2 <- ggdraw(p.u + theme(legend.position="none") + labs(title = title[1])) +
          draw_label(bquote(paste(p == .(round(u.adonis[row,"Pr(>F)"],3)),"; ", R^2 == .(round(u.adonis[row,"R2"],3)))),
                     x = Inf, y = Inf, hjust = 1.05, vjust = 0.3, size = 14) +
          draw_label(bquote(paste(p == .(round(u.adonis[2,"Pr(>F)"],3)),"; ", R^2 == .(round(u.adonis[2,"R2"],3)))),
                     x = Inf, y = Inf, hjust = 1.05, vjust = 1.15, size = 14)
        
        p.w2 <- ggdraw(p.w + theme(legend.position="none") + labs(title = title[1])) +
          draw_label(bquote(paste(p == .(round(w.adonis[row,"Pr(>F)"],3)),"; ", R^2 == .(round(w.adonis[row,"R2"],3)))),
                     x = Inf, y = Inf, hjust = 1.05, vjust = 0.3, size = 14) +
          draw_label(bquote(paste(p == .(round(w.adonis[2,"Pr(>F)"],3)),"; ", R^2 == .(round(w.adonis[2,"R2"],3)))),
                     x = Inf, y = Inf, hjust = 1.05, vjust = 1.15, size = 14)
      }
      grouped_plots <- plot_grid(p.u2, p.w2, ncol = 2, align = "hv")
      
      p <- plot_grid(arrangeGrob(nullGrob(), grouped_plots, nullGrob(), ncol=1, heights=c(.1,2,.2)),
                     legend, ncol=2, rel_widths = c(.3, .05))
      
    } else {
      p <- ggdraw(p) + labs(title = title[1]) +
        draw_label(bquote(paste(p == .(round(u.adonis[row,"Pr(>F)"],3)),"; ", R^2 == .(round(u.adonis[row,"R2"],3)))),
                   x = Inf, y = Inf, hjust = 1.1, vjust = 1, size = 14)
      
      if(!is.null(var2) && inter == FALSE){
        p <- p +
          draw_label(bquote(paste(p == .(round(u.adonis[2,"Pr(>F)"],3)),"; ", R^2 == .(round(u.adonis[2,"R2"],3)))),
                            x = Inf, y = Inf, hjust = 1.1, vjust = 1.8, size = 14)
      }
    }
  } else {beta_stats <- list()}
  if (var1 == "Sample") {
    p <- p + theme(legend.position = "none")
    p <- ggplotly(p, tooltip = c("Sample", "Axis.1", "Axis.2"))
    beta_stats <- p
  } else {beta_stats[["plot"]] <- p}
  if (isTRUE(print)) {
    print(p)
  }
  return(beta_stats)
}
