#' Calculate and plot alpha diversity results
#'
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr group_by_at
#' @importFrom dplyr summarize
#' @importFrom dplyr select_at
#' @importFrom dplyr arrange
#' @importFrom dplyr n
#' @importFrom magrittr %>%
#' @importFrom ggpubr ggboxplot
#' @importFrom ggpubr stat_compare_means
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom ggpubr rotate_x_text
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 expansion
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 margin
#' @importFrom vegan specnumber
#' @importFrom vegan diversity
#' @importFrom rstatix kruskal_test
#' @importFrom rstatix wilcox_test
#' @importFrom rstatix add_significance
#' @importFrom rstatix add_xy_position
#' @importFrom QsRutils make_letter_assignments
#' @param data (Required). The list object imported with \code{\link{import_files}}.
#' @param var1 (Required). Character string. variable to evaluate.
#' @param var2 (Optional). Character string.
#'  Secondary variable only for alpha diversity stats.
#' @param div (Required). Logical. Default \code{TRUE}.
#'  Include diversity analysis with the index passed by the \code{idx} parameter.
#' @param idx (Required). Character string. Default \code{Shannon}.
#'  If \code{div} is \code{TRUE}, it specified the diversity index to be used along with richness.
#' @param facet.var2 (Required). Logical. Default \code{FALSE}.
#'  Pass \code{var2} to \code{\link{facet_wrap}} in the box-plots.
#' @param stats (Required). Logical. Default \code{TRUE}.
#'  Include overall Kruskal-Wallis results in the plot.
#' @param pairwise (Required). Logical. Default \code{TRUE}.
#'  Include overall pairwise Wilcoxon results in the plot.
#'  If \code{facet.var2} is \code{TRUE}, pairwise results will be represented by brackets.
#' @param levels.1 (Optional). Character vector. Sort var1.
#' @param levels.2 (Optional). Character vector. Sort var2.
#' @param xlab (Optional). Character string. Title for axis x in the plot.
#' @param print (Required). Logical. Default \code{TRUE}. Print the plot.
#' @return Returns a list object and prints a plot.
#' @usage plot_alpha(data=NULL, var1=NULL, var2=NULL, div=TRUE, idx="Shannon",
#'  facet.var2=FALSE, stats=TRUE, pairwise=TRUE, levels.1=NULL, levels.2=NULL, xlab=NULL, print=TRUE)
#' @examples
#' library(freyR)
#' # Import QIIME2 artifacts and metadata
#' files <- import_files(metadata="Path/metadata.tsv", table_norm="Path/table_norm.qza", taxonomy="Path/taxonomy.qza")
#'
#' # Calculate and plot alpha diversity results
#' alpha <- plot_alpha(files, var1="Treatment")
#' @export
plot_alpha <- function(data, var1, var2=NULL, div=TRUE, idx="Shannon", facet.var2=FALSE,
                       stats=TRUE, pairwise=TRUE, levels.1=NULL, levels.2=NULL, xlab=NULL, print=TRUE) {
  #Subset data components
  metadata <- data[["metadata"]]
  table_norm <- data[["table_norm"]]

  #Calculate alpha diversity
  richness <- data.frame(Score = specnumber(t(table_norm)), Index = "Richness", SampleID = colnames(table_norm))
  alpha <- merge(richness, metadata, by = "SampleID", all.y = F)
  if (div == TRUE) {
    diversity <- data.frame(Score = diversity(t(table_norm), index = tolower(idx), base = 2),
                            Index = idx, SampleID = colnames(table_norm))
    alpha <- merge(rbind(richness, diversity), metadata, by = "SampleID", all.y = F) %>% arrange(Index)
  }
  names(alpha)[1] <- "Sample"
  alpha[alpha == ''] <- NA
  alpha[,colnames(alpha) == var1] <- factor(alpha[,colnames(alpha) == var1])

  #Stat tables
  stats_table <- function(x, var){
    x %>%
      filter(!is.na(.data[[var]])) %>%
      group_by_at(c(var, "Index")) %>%
      summarize(
        N = n(),
        mean = round(mean(Score),2),
        median = round(median(Score),2),
        min = round(min(Score),2),
        max = round(max(Score),2),
        lower.q = round(quantile(Score, 0.25),2),
        upper.q = round(quantile(Score, 0.75),2)) %>%
      arrange(Index)
  }

  alpha_table_v1 <- stats_table(alpha, var1)

  alpha_mean <- alpha %>% group_by(Index) %>% summarize(mean = mean(Score))

  if (var1 == "Sample") {
    stats = FALSE
    alpha_stats = list(var1 = alpha_table_v1)
  } else if (div == TRUE) {
    alpha_table <- data.frame(Sample = alpha$Sample, Richness = alpha$Score[alpha$Index == "Richness"],
                              Shannon = alpha$Score[alpha$Index == "Shannon"]) %>% distinct()
    alpha_stats = list(var1 = alpha_table_v1, Sample = alpha_table)
  } else {
    alpha_table <- data.frame(Sample = alpha$Sample, Richness = alpha$Score[alpha$Index == "Richness"]) %>% distinct()
    alpha_stats = list(var1 = alpha_table_v1, Sample = alpha_table)
  }

  if (!is.null(var2)) {
    alpha[,colnames(alpha) == var2] <- factor(alpha[,colnames(alpha) == var2])
    color = var2
    alpha_table_v2 <- stats_table(alpha, var2)
    alpha_stats <- list(var1 = alpha_table_v1, var2= alpha_table_v2, Sample = alpha_table)
  } else {color = var1}

  #Base graph
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                 "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

  if (!is.null(levels.1)) {
    alpha[,colnames(alpha) == var1] <- factor(alpha[,colnames(alpha) == var1], levels = levels.1)
  }
  if (!is.null(levels.2)) {
    alpha[,colnames(alpha) == var2] <- factor(alpha[,colnames(alpha) == var2], levels = levels.2)
  }
  alpha_noNA <- alpha %>% filter(!is.na(.data[[var1]]))

  if (nrow(unique(alpha_noNA[var1])) < 5 || (div == TRUE && facet.var2 == TRUE)) {ncol = 2} else {ncol = 1}

  p <- ggboxplot(alpha_noNA, x = var1, y = "Score", color = color, add = "jitter") +
    facet_wrap("Index", scales = "free_y", ncol = ncol) + labs(x = xlab) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    geom_hline(data= alpha_mean, aes(yintercept=mean), linetype = 2) + # Add horizontal line at base mean
    theme_light()

  p <- p + theme(axis.title.y = element_blank(), text=element_text(size=12), plot.margin = margin(1,1,1,1.5, "cm")) +
    rotate_x_text(angle = 45) +
    scale_colour_manual(values=rep(cbPalette, times = 25))

  if (is.null(var2) | facet.var2 == TRUE) {
    p <- p + facet_wrap(c("Index", var2), scales = "free_y", ncol = ncol) + theme(legend.position = "none")
  }

  #Add stats to the plot
  if (stats == FALSE) {
    alpha_stats[["plot"]] <- p
    print(p)
    return(alpha_stats)
  } else {
    options(scipen = 999)
    stat <- alpha_noNA %>%
      group_by(Index) %>%
      kruskal_test(as.formula(paste("Score ~", var1)))
    alpha_stats[["stat"]] <- stat

    if (length(unique(alpha_noNA[,colnames(alpha_noNA)==var1]))>2 & pairwise == TRUE){
      stat.pairwise <- alpha_noNA %>%
        group_by(Index) %>%
        wilcox_test(as.formula(paste("Score ~", var1)), p.adjust.method = "BH") %>%
        add_significance() %>%
        add_xy_position(x = var1, scales = "free", fun = "max")

      alpha_stats[["stat.pairwise.r"]] <- filter(stat.pairwise, stat.pairwise$Index == "Richness")
      if (length(unique(alpha_noNA[,colnames(alpha_noNA)==var1]))<5||facet.var2 == TRUE) {
        #Alpha-plot with pairs indicated using brackets
        if(!all(stat.pairwise$p.adj.signif[stat.pairwise$Index == "Richness"] == "ns")) {
          p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == "Richness",], method = "kruskal.test",
                                      label.y = (max(stat.pairwise$y.position)*1.1), label.x = 1.3) +
            stat_pvalue_manual(stat.pairwise, hide.ns = TRUE, xmin = "group1",
                               step.increase = 0.02, step.group.by = "Index", label = "p.adj.signif")
        } else {
          p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == "Richness",], method = "kruskal.test",
                                      label.y = max(alpha_noNA$Score)*1.12, label.x = 1.3)
        }
      } else {
        #Alpha-plot with pairs indicated using compact letter display
        alpha.max <- alpha_noNA %>% group_by_at(c("Index", var1)) %>% summarize(max = max(Score), min = min(Score))

        ptt.rslt.r <- with(alpha_noNA[alpha_noNA$Index=="Richness",],
                           pairwise.wilcox.test(Score, alpha_noNA[alpha_noNA$Index=="Richness",colnames(alpha_noNA)==var1],
                                                p.adjust.method = "BH"))

        ltrs.r <- make_letter_assignments(ptt.rslt.r)

        ltr_df <- data.frame(x = rep(c(1:length(ltrs.r$Letters)),2), alpha.max, ltrs = ltrs.r$Letters)

        p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == "Richness",], method = "kruskal.test",
                                    label.y = (max(ltr_df$max[ltr_df$Index == "Richness"])*1.3), label.x = 1.3)
        if (div == FALSE) {
          p <- p + geom_text(data = ltr_df, aes(x=x, y=max, label=ltrs), nudge_y = median(ltr_df$max)*0.1, fontface = "bold")}
      }
    } else if (is.null(var2) | facet.var2 == TRUE) {
      p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == "Richness",], method = "kruskal.test",
                                  label.x = 1.3, label.y = max(alpha_noNA$Score)*1.03)
    } else {
      p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == "Richness",], method = "kruskal.test",
                                  label.x = 1.3, label.y = max(alpha_noNA$Score)*1.12) +
        stat_compare_means(aes(group = .data[[var2]]), label = "p.signif")
    }
    #Add a second plot for the selected diversity index
    if (div == TRUE){
      if (length(unique(alpha_noNA[,colnames(alpha_noNA)==var1]))>2 & pairwise == TRUE){
        alpha_stats[["stat.pairwise.s"]] <- filter(stat.pairwise, stat.pairwise$Index == idx)
        if (length(unique(alpha_noNA[,colnames(alpha_noNA)==var1]))<5 || facet.var2 == TRUE) {
          #Alpha-plot with pairs indicated using brackets
          if(!all(stat.pairwise$p.adj.signif == "ns")) {
            p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == idx,], method = "kruskal.test",
                                        label.y = (max(stat.pairwise[stat.pairwise$Index == idx,]$y.position)*1.08), label.x = 1.3)
          } else {
            p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == idx,], method = "kruskal.test",
                                        label.y = max(alpha_noNA[alpha_noNA$Index == idx,]$Score)*1.12, label.x = 1.3)
          }
        } else {
          #Alpha-plot with pairs indicated using compact letter display
          ptt.rslt.s <- with(alpha_noNA[alpha_noNA$Index==idx,],
                             pairwise.wilcox.test(Score, alpha_noNA[alpha_noNA$Index==idx,colnames(alpha_noNA)==var1],
                                                  p.adjust.method = "BH"))

          ltrs.s <- make_letter_assignments(ptt.rslt.s)

          ltr_df <- data.frame(x = rep(c(1:length(ltrs.r$Letters)),2), alpha.max, ltrs = c(ltrs.r$Letters, ltrs.s$Letters)) %>%
            group_by(Index) %>%
            mutate(nudge_y = median(ifelse(max != min, max-min, max))*0.2)

          p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == idx,], method = "kruskal.test",
                                      label.y = (max(ltr_df$max[ltr_df$Index == idx])*1.5), label.x = 1.3) +
            geom_text(data = ltr_df, aes(x=x, y=max, label=ltrs), nudge_y = ltr_df$nudge_y, fontface = "bold")
        }
      } else {
        p <- p + stat_compare_means(data = alpha_noNA[alpha_noNA$Index == idx,], method = "kruskal.test",
                                    label.x = 1.3, label.y = max(alpha_noNA[alpha_noNA$Index == idx,]$Score)*1.12)
      }
    }
    alpha_stats[["plot"]] <- p
    print(p)
    return(alpha_stats)
  }
}
