######################################################################################################################
# R script to assign cell type to cells in scRNA-seq data generated on 10xGenomics Chromium platform                 #
# This version uses Haemopedia, with classifications revised April 2019                                              #
# Creates a subdirectory of the input directory                                                                      #
# called CellMatch.Haemopedia.April2019, containing the following:                                                   #
# 1. CellMatch.Haemopedia.lineage.%s.tsv                                                                             #
# 2. CellMatch.Haemopedia.lineage.partek.%s.tsv                                                                      #
# 3. CellMatch.Haemopedia.lineage.loupe.%s.tsv                                                                       #
# Requires Seurat V3, R >= 3.5                                                                                       #
# Developed by: Allegra Petti (allegra.petti@wustl.edu)                                                              #
# Reformatted by: Sai Mukund (saimukund@wustl.edu)                                                                   #
# The code should be run as follows: CellMatch.Haemopedia.190507.r Samplename /path/to/input/dir /path/to/output/dir #
######################################################################################################################

#Checking availability of packages
if(!require("foreach"))
{
   print ("\n\nPackage foreach doesn't exist, so I'm installing it!\n\n")
   install.packages("foreach")
   library(foreach)
}
if(!require("doParallel"))
{
   print ("\n\nPackage doParallel doesn't exist, so I'm installing it!\n\n")
   install.packages("doParallel")
   library(doParallel)
}

library("Seurat");
library("dplyr");
library(foreach);
library(doParallel);
registerDoParallel(cores=12);

# read from command line (see examples below):
args = commandArgs(trailingOnly=TRUE);
if(length(args) != 6)
{
   stop("\n\nAll the arguments were not provided.\n\nUser must provide sample name, path to input root directory and path to output directory.\n\nThe code should be run as follows: CellMatch_Haemopedia.r Samplename /path/to/input/dir /path/to/output/dir\n\n")
}

upn = args[1]; # upn or sample name
matrix.dir = args[2]; # root directory for input
output.main = args[3]; # path for outputs
signature.file = args[4]; #path for the lineage dataset
mincell = args[5]
minfeat = args[6]


dir.create(file.path(output.main))
date = gsub("2019-","19",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

# load the data
print ("Loading data...")
data.10x <- Read10X(matrix.dir);

# Initialize the Seurat object with the raw (non-normalized data)
print ("Creating Seurat object...")

if(mincell>=3 && minfeat>=10){
 scrna <- CreateSeuratObject(counts=data.10x, min.cells=mincell, min.features=minfeat, project=sprintf("%s.%s",upn,date))
} else{
 scrna <- CreateSeuratObject(counts=data.10x, project=sprintf("%s.%s",upn,date))
}

# Normalize the data
print ("Normalizing data...")
scrna <- NormalizeData(object = scrna, normalization.method = "LogNormalize", scale.factor = 1e4)
genes <- rownames(x = scrna); # genes in normalized seurat object
cells <- Cells(scrna); # cells in normalized seurat object

# load the haemopedia gene signatures:
lineage.header = scan(signature.file, sep='\t', nlines=1, what=character(), fill=TRUE);
lineage.list=lineage.header[4:length(lineage.header)]
haem.raw = read.table(signature.file, sep='\t', skip=1, header=FALSE, fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE);
names(haem.raw) = lineage.header;
haem.unsorted = haem.raw[which(haem.raw$Gene %in% genes),]; # store the haemopedia data for genes that are also in the seurat object

### extract haemopedia gene data from normalized seurat object:
g.ind = which(genes %in% haem.unsorted$Gene);
scrna.subset = FetchData(object=scrna, vars=genes[g.ind]); # now rows are cells, and columns are genes
scrna.subset = t(scrna.subset); # Transpose so rows are genes, columns are cells. this is scRNAseq data

# sort the haemopedia data so the genes are in the same order as in scrna.subset, above
new.indices = match(rownames(scrna.subset),haem.unsorted$Gene);
haem.data <- haem.unsorted[new.indices,]; # this is haemopedia data used below

# set up the data frame to hold the results:
results <- foreach(i=1:length(cells), .combine=rbind) %dopar% { # for each cell
#results <- foreach(i=1:10, .combine=rbind) %dopar% { # for each cell
	 if (i%%100 == 0) {
    print(i);
  }
  profile = as.vector(scrna.subset[,which(colnames(scrna.subset)==cells[i])]);
corrs = vector(length=ncol(haem.data)-3, mode="numeric"); # list of correlations, one for each signature
   k = 0;
   for (j in 4:ncol(haem.data)) { # for each column/lineage in the signature matrix
     k = k+1;
     corrs[k] = cor(profile,haem.data[,j],method="spearman");
   }
   order.ind = order(corrs, decreasing=TRUE); # indices of the corrs vector in decreasing order
# get indices and values for the top three predictions:
   i1 = order.ind[1];
   i2 = order.ind[2];
   i3 = order.ind[3];
   max1 = corrs[i1];
   max2 = corrs[i2];
   max3 = corrs[i3];
   lineage1 = as.character(lineage.list[i1]);
   lineage2 = as.character(lineage.list[i2]);
   lineage3 = as.character(lineage.list[i3]);
   linlist = c(lineage1, lineage2, lineage3);
   if (length(unique(linlist)) == 3) {
      LINEAGE = "NP"; # no prediction possible
   } else {
      lintab <- table(linlist); # make a count table
      maxlin <- max(lintab); # max count in the table
      LINEAGE <- names(lintab[which(lintab == maxlin)]); # name of lineage w/ max count
   }
#   print(LINEAGE)
   # put the results into a data table:
   results.temp = data.frame(as.character(cells[i]), LINEAGE, lineage1, lineage2, lineage3, max1, max2, max3);
} # end for each cell

colnames(results) <- c("CellBarcode", "LINEAGE", "lineage1", "lineage2", "lineage3", "max1", "max2", "max3");
write.table(results, file=paste0(output.main,"/",upn,"_CellMatch.Haemopedia.lineage.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# loupe browser output
loupe <- results[,c(1,2)];
colnames(loupe) <- c("Barcode","Lineage");
if (length(grep("-\\d",loupe$Barcode)) == 0) { # if the cells do not have flags on them, add a one-sample flag
   loupe$Barcode <- sprintf("%s-1",loupe$Barcode);
}
write.table(loupe, file=paste0(output.main,"/",upn,"_CellMatch.Haemopedia.lineage_loupe.tsv"), sep=",", quote=FALSE, row.names=FALSE)

# now rename for Partek Output:
colnames(loupe) <- c("Cell name","Cell Type");
write.table(loupe, file=paste0(output.main,"/",upn,"_CellMatch.Haemopedia.lineage_partek.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
