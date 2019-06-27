#base packages
install.packages(c("tidyverse","ggplot2","Hmisc","plotrix","png","gplots","dplyr","Seurat","foreach","doParallel"),repo=paste0("https://mran.microsoft.com/snapshot/",format(Sys.Date(), format="%Y-%m-%d")))

#bioconductor packages
#this is modified from from https://raw.githubusercontent.com/Bioconductor/bioc_docker/master/src/core/install.R.in
pkgs <- c("bsseq","data.table","methylKit","readr","sqldf","stringr", "scater", "slingshot")

ap.db <- available.packages(contrib.url(BiocManager::repositories()))
ap <- rownames(ap.db)

pkgs_to_install <- pkgs[pkgs %in% ap]

BiocManager::install(pkgs_to_install, update=FALSE, ask=FALSE)

# just in case there were warnings, we want to see them
# without having to scroll up:
warnings()

if (!is.null(warnings()))
{
    w <- capture.output(warnings())
    if (length(grep("is not available|had non-zero exit status", w)))
        quit("no", 1L)
}

suppressWarnings(BiocManager::install(update=TRUE, ask=FALSE))
