# Installation
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
devtools::install_git("https://gitlab.fera.co.uk/mfernand/freyR")