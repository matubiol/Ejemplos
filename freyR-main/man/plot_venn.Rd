% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_venn.R
\name{plot_venn}
\alias{plot_venn}
\title{Generates Venn diagrams and provides list of taxa features}
\usage{
plot_venn(data=NULL, norm=FALSE, taxrank=c("Species","Genus",...), var1=NULL, var2=NULL, freq=NULL,
 comb=FALSE, upset=FALSE, names=NULL, title=NULL, write=FALSE, print=TRUE)
}
\arguments{
\item{data}{(Required). The list object imported with \code{\link{import_files}}.}

\item{norm}{(Required). Logical.
If \code{TRUE}, normalized table is used.}

\item{taxrank}{(Optional). Character string. Taxonomic level to represent in the diagram.}

\item{var1}{(Required). Character string. Variable that determines the initial level of sample subsets.
If \code{var2} is \code{NULL}, the number of circles in the Venn diagram will match the number of groups in \code{var1}.}

\item{var2}{(Optional). Character string. Variable that determines the second level of sample subsets.}

\item{freq}{(Optional). Numeric.
The percentage of the number of samples for each \code{var1/var2} in which a 'feature' must appear.}

\item{comb}{(Required). Logical. Default \code{FALSE}.
Combine \code{var1} and \code{var2} in diagram. If the number of groups > 4, an \code{\link{upset}} plot will be created.}

\item{NArm}{(Optional). Logical. Default \code{TRUE}. CAUTION.
The decision to prune (or not) taxa for which you lack categorical data could have a large effect on downstream analysis.}

\item{upset}{(Required). Logical. Default \code{FALSE}. Force an \code{\link{upset}} plot, even if the number of groups < 5.}

\item{names}{(Optional). Character string. This allows customization of the group names in the diagram.}

\item{sort}{(Optional). Character vector. Sort groups. If \code{names} is not \code{NULL}, use new names.}

\item{title}{(Optional). Character string. This add a title to the diagram.}

\item{write}{(Required). Logical. Default \code{FALSE}.}

\item{print}{(Required). Logical. Default \code{TRUE}. Print the plot.
If \code{TRUE}, the list of features within each area of the diagram is written to an Excel file in the working directory.}
}
\value{
Returns a list object and prints a diagram.
}
\description{
Generates Venn diagrams and provides list of taxa features
}
\examples{
library(freyR)
# Import QIIME2 artifacts and metadata
files <- import_files(metadata="Path/metadata.tsv", table="Path/table.qza", taxonomy="Path/taxonomy.qza")

# Plot venn and get the features of each group.
venn <- plot_venn(files, taxrank = "Species", var1 = "Treatment", var2 = "Time", freq = 0.3, comb = FALSE)
}
