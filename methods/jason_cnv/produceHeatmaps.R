# Options
myDir           <- "/Users/jason/Desktop/OneDrive - University of Minnesota/BICB/Thesis/CNV/"
ensemblGenesDir <- "/Users/jason/Desktop/OneDrive - University of Minnesota/BICB/Thesis/Data/Human/ensembl_genes/"
patientDir      <- paste("/Users/jason/Desktop/OneDrive - University of Minnesota/BICB/Thesis/Data/Patients/",sep="")
patient         <- "p59"
# patients <- c("p59","p60","p62","p64","p67","p69","p71","p72","p75","p76", "p77", "p78", "p82", "p83", "p85")

# initialize
setwd(myDir)
library(ggplot2)
library(data.table)
library(readxl)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(gplots)

# Ingest & annotate geneData from human reference data
chr1 <- read_excel(paste(ensemblGenesDir,"Chromosome 1.xls",sep = ""))
chr2 <- read_excel(paste(ensemblGenesDir,"Chromosome 2.xls",sep = ""))
chr3 <- read_excel(paste(ensemblGenesDir,"Chromosome 3.xls",sep = ""))
chr4 <- read_excel(paste(ensemblGenesDir,"Chromosome 4.xls",sep = ""))
chr5 <- read_excel(paste(ensemblGenesDir,"Chromosome 5.xls",sep = ""))
chr6 <- read_excel(paste(ensemblGenesDir,"Chromosome 6.xls",sep = ""))
chr7 <- read_excel(paste(ensemblGenesDir,"Chromosome 7.xls",sep = ""))
chr8 <- read_excel(paste(ensemblGenesDir,"Chromosome 8.xls",sep = ""))
chr9 <- read_excel(paste(ensemblGenesDir,"Chromosome 9.xls",sep = ""))
chr10 <- read_excel(paste(ensemblGenesDir,"Chromosome 10.xls",sep = ""))
chr11 <- read_excel(paste(ensemblGenesDir,"Chromosome 11.xls",sep = ""))
chr12 <- read_excel(paste(ensemblGenesDir,"Chromosome 12.xls",sep = ""))
chr13 <- read_excel(paste(ensemblGenesDir,"Chromosome 13.xls",sep = ""))
chr14 <- read_excel(paste(ensemblGenesDir,"Chromosome 14.xls",sep = ""))
chr15 <- read_excel(paste(ensemblGenesDir,"Chromosome 15.xls",sep = ""))
chr16 <- read_excel(paste(ensemblGenesDir,"Chromosome 16.xls",sep = ""))
chr17 <- read_excel(paste(ensemblGenesDir,"Chromosome 17.xls",sep = ""))
chr18 <- read_excel(paste(ensemblGenesDir,"Chromosome 18.xls",sep = ""))
chr19 <- read_excel(paste(ensemblGenesDir,"Chromosome 19.xls",sep = ""))
chr20 <- read_excel(paste(ensemblGenesDir,"Chromosome 20.xls",sep = ""))
chr21 <- read_excel(paste(ensemblGenesDir,"Chromosome 21.xls",sep = ""))
chr22 <- read_excel(paste(ensemblGenesDir,"Chromosome 22.xls",sep = ""))
chrX <- read_excel(paste(ensemblGenesDir,"Chromosome X.xls",sep = ""))

geneData <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX)

colnames(geneData) <- c("ENSEMBL", "chr", "start", "end", "geneName")

rm(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX)

geneData$arm <- "NA"
geneData$arm[geneData$chr == 1 & geneData$start < 125000000] <- 'p'
geneData$arm[geneData$chr == 1 & geneData$start > 125000000] <- 'q'
geneData$arm[geneData$chr == 2 & geneData$start < 93300000] <- 'p'
geneData$arm[geneData$chr == 2 & geneData$start > 93300000] <- 'q'
geneData$arm[geneData$chr == 3 & geneData$start < 91000000] <- 'p'
geneData$arm[geneData$chr == 3 & geneData$start > 91000000] <- 'q'
geneData$arm[geneData$chr == 4 & geneData$start < 50400000] <- 'p'
geneData$arm[geneData$chr == 4 & geneData$start > 50400000] <- 'q'
geneData$arm[geneData$chr == 5 & geneData$start < 48400000] <- 'p'
geneData$arm[geneData$chr == 5 & geneData$start > 48400000] <- 'q'
geneData$arm[geneData$chr == 6 & geneData$start < 61000000] <- 'p'
geneData$arm[geneData$chr == 6 & geneData$start > 61000000] <- 'q'
geneData$arm[geneData$chr == 7 & geneData$start < 59900000] <- 'p'
geneData$arm[geneData$chr == 7 & geneData$start > 59900000] <- 'q'
geneData$arm[geneData$chr == 8 & geneData$start < 45600000] <- 'p'
geneData$arm[geneData$chr == 8 & geneData$start > 45600000] <- 'q'
geneData$arm[geneData$chr == 9 & geneData$start < 49000000] <- 'p'
geneData$arm[geneData$chr == 9 & geneData$start > 49000000] <- 'q'
geneData$arm[geneData$chr == 10 & geneData$start < 40200000] <- 'p'
geneData$arm[geneData$chr == 10 & geneData$start > 40200000] <- 'q'
geneData$arm[geneData$chr == 11 & geneData$start < 53700000] <- 'p'
geneData$arm[geneData$chr == 11 & geneData$start > 53700000] <- 'q'
geneData$arm[geneData$chr == 12 & geneData$start < 35800000] <- 'p'
geneData$arm[geneData$chr == 12 & geneData$start > 35800000] <- 'q'
geneData$arm[geneData$chr == 13 & geneData$start < 17900000] <- 'p'
geneData$arm[geneData$chr == 13 & geneData$start > 17900000] <- 'q'
geneData$arm[geneData$chr == 14 & geneData$start < 17600000] <- 'p'
geneData$arm[geneData$chr == 14 & geneData$start > 17600000] <- 'q'
geneData$arm[geneData$chr == 15 & geneData$start < 19000000] <- 'p'
geneData$arm[geneData$chr == 15 & geneData$start > 19000000] <- 'q'
geneData$arm[geneData$chr == 16 & geneData$start < 36600000] <- 'p'
geneData$arm[geneData$chr == 16 & geneData$start > 36600000] <- 'q'
geneData$arm[geneData$chr == 17 & geneData$start < 24000000] <- 'p'
geneData$arm[geneData$chr == 17 & geneData$start > 24000000] <- 'q'
geneData$arm[geneData$chr == 18 & geneData$start < 17200000] <- 'p'
geneData$arm[geneData$chr == 18 & geneData$start > 17200000] <- 'q'
geneData$arm[geneData$chr == 19 & geneData$start < 26500000] <- 'p'
geneData$arm[geneData$chr == 19 & geneData$start > 26500000] <- 'q'
geneData$arm[geneData$chr == 20 & geneData$start < 27500000] <- 'p'
geneData$arm[geneData$chr == 20 & geneData$start > 27500000] <- 'q'
geneData$arm[geneData$chr == 21 & geneData$start < 13200000] <- 'p'
geneData$arm[geneData$chr == 21 & geneData$start > 13200000] <- 'q'
geneData$arm[geneData$chr == 22 & geneData$start < 14700000] <- 'p'
geneData$arm[geneData$chr == 22 & geneData$start > 14700000] <- 'q'
geneData$arm[geneData$chr == "X" & geneData$start < 60600000] <- 'p'
geneData$arm[geneData$chr == "X" & geneData$start > 60600000] <- 'q'

# read in files and re-label headers
barcodes.tsv <- read.csv(paste(patientDir,patient,"/barcodes.tsv", sep=""), header = FALSE)
colnames(barcodes.tsv) <- "barcodes"

genes.tsv <- read.csv(paste(patientDir,patient,"/genes.tsv", sep=""), header = FALSE, sep = "\t")
colnames(genes.tsv) <- c("ENSEMBL","Names")

matrix.mtx <- read.csv(paste(patientDir,patient,"/matrix.mtx", sep=""), header = FALSE, skip = 3, sep = " ")
colnames(matrix.mtx) <- c("geneIndex", "cellIndex", "UMI")

if ( file.exists(paste(patientDir,patient,"/cluster_annotation.xlsx", sep="")) ) 
{ 
  clusterDataExists <- TRUE 
} else { 
  clusterDataExists <- FALSE 
}

if ( clusterDataExists )
{ 

#  cluster_data.csv <- read.csv(paste(patientDir,patient,"/cluster_data.csv", sep=""), header = TRUE, skip = 0, sep = ",") 
  cluster_data.csv <- read_excel(paste(patientDir,patient,"/cluster_annotation.xlsx", sep=""))  
  
  colnames(cluster_data.csv) <- c("origIndex", "filterIndex", "barcode", 
                                  "clustNum", "clustColor", 
                                  "majorCellType", "majorCellTypeColor", 
                                  "minorCellType", "minorCellTypeColor")
} else { 
  annot1 <- read.csv(paste(patientDir,patient,"/annot-01.csv", sep=""), header = TRUE, sep = ",")
  annot2 <- read.csv(paste(patientDir,patient,"/annot-02.csv", sep=""), header = TRUE, sep = ",") 
  colnames(annot1) <- c("barcode", "clustNum")
  colnames(annot2) <- c("barcode", "cellType")
}

# add ensemble index from genes.tsv file
geneData$ensemblIndex <- match(geneData$ENSEMBL, genes.tsv$ENSEMBL)

# use indices of genes to match to UMI values in matrix
masterMatrix <- merge(geneData, matrix.mtx, by.x="ensemblIndex", by.y="geneIndex", all.x=FALSE, all.y = FALSE)


# REDUCE TO HOUSEKEEPING ONLY
housekeepingGenes <- read.csv(file = "/Users/jason/Desktop/OneDrive - University of Minnesota/BICB/Thesis/Data/Human/3804_housekeeping_gene_list.txt", sep = "\t", header = FALSE)
colnames(housekeepingGenes) <- c("gene", "otherName")
housekeepingGenes <- trimws(as.character(housekeepingGenes$gene))
masterMatrix <- masterMatrix[masterMatrix$geneName %in% housekeepingGenes,]


# add barcodes, cluster numbers, and cell types to masterMatrix
masterMatrix$barcode <- barcodes.tsv[masterMatrix$cellIndex,]
if ( clusterDataExists )
{
  masterMatrix$clustNum <- cluster_data.csv$clustNum[match(masterMatrix$barcode, cluster_data.csv$barcode)] 
  masterMatrix$cellType <- cluster_data.csv$minorCellType[match(masterMatrix$barcode, cluster_data.csv$barcode)]
} else {
  masterMatrix$clustNum <- annot1$clustNum[match(masterMatrix$barcode, annot1$barcode)]
  masterMatrix$cellType <- annot2$cellType[match(masterMatrix$barcode, annot2$barcode)]
}
masterMatrix$chrArm <- paste(masterMatrix$chr, masterMatrix$arm, sep = "")

### Turn UMI counts into fractions of total cell UMI by dividing each UMI by sum of all UMIs for that cell
totCellUMI <- aggregate( UMI ~ barcode, data = masterMatrix, FUN = sum )
colnames(totCellUMI) <- c("barcode", "sumUMI")
masterMatrix <- merge(masterMatrix, totCellUMI, by.x="barcode", by.y="barcode", all.x=FALSE, all.y = FALSE)
rm(totCellUMI)
masterMatrix$UMIfrac <- masterMatrix$UMI / masterMatrix$sumUMI

# sum value per chrArm  
chrArmFracSums <- aggregate( UMIfrac ~ barcode + chrArm, data = masterMatrix, FUN = sum  )
# add cell type data
chrArmFracSums <- merge(x = chrArmFracSums, y = cluster_data.csv, by = "barcode", all.x = FALSE, all.y=FALSE)

#### calculate normalized values
cellCenter  <- aggregate( formula = UMIfrac ~ chrArm, data = chrArmFracSums[chrArmFracSums$minorCellType != "epithelial" & chrArmFracSums$minorCellType != "epithelial_cycling",], FUN = mean )
cellScale   <- aggregate( formula = UMIfrac ~ chrArm, data = chrArmFracSums[chrArmFracSums$minorCellType != "epithelial" & chrArmFracSums$minorCellType != "epithelial_cycling",], FUN = sd )
#cellCenter <- aggregate( formula = UMIfrac ~ chrArm, data = chrArmFracSums[chrArmFracSums$clustNum %in% c(7),],            FUN = mean )
#cellScale  <- aggregate( formula = UMIfrac ~ chrArm, data = chrArmFracSums[chrArmFracSums$clustNum %in% c(7),],            FUN = sd )

colnames(cellCenter)   <- c("chrArm", "meanUMI")
colnames(cellScale)    <- c("chrArm", "sdUMI")

chrArmFracSums <- merge(x = chrArmFracSums, y = cellCenter,   by = "chrArm", all.x = FALSE, all.y = FALSE)
chrArmFracSums <- merge(x = chrArmFracSums, y = cellScale,    by = "chrArm", all.x = FALSE, all.y = FALSE)

chrArmFracSums$UMInorm   <- ( chrArmFracSums$UMIfrac - chrArmFracSums$meanUMI   ) / chrArmFracSums$sdUMI


### create heatmap normalized on all non-epithelial cells
######################
info <- paste(patient, "total fraction of UMI counts per cell, centered/scaled by chrArm on all non epithelial cells", sep = "; ")

# convert data to wide matrix
heatmapMatrix <- chrArmFracSums[,c("chrArm", "barcode", "UMInorm")] 
heatmapMatrix <- dcast(heatmapMatrix, barcode~chrArm, value.var="UMInorm") 

# add cell type, clust num, and color data
heatmapMatrix$cellType <- cluster_data.csv$minorCellType[match(heatmapMatrix$barcode, cluster_data.csv$barcode, nomatch = "NA")]
heatmapMatrix$clustNum <- cluster_data.csv$clustNum[match(heatmapMatrix$barcode, cluster_data.csv$barcode, nomatch = "NA")]
heatmapMatrix$color    <- cluster_data.csv$minorCellTypeColor[match(heatmapMatrix$barcode, cluster_data.csv$barcode, nomatch = "NA")]

# fix rownames
rownames(heatmapMatrix) <- heatmapMatrix[,"barcode"] 

# order y-axis (cells) by cellType, then clustNum
heatmapMatrix <- heatmapMatrix[order(heatmapMatrix$cellType, heatmapMatrix$clustNum),]

# create row seperators based on cluster numbers
rowSeps <- c(1,1+which(diff(heatmapMatrix$clustNum)!=0))

# save heatmap colors
heatmapColors <- heatmapMatrix$color
heatmapColors[heatmapColors==""] <- NA

# make group row labels
groupLabels <- rep("", nrow(heatmapMatrix))
groupLabels[rowSeps] <- paste(heatmapMatrix$clustNum[rowSeps], heatmapMatrix$cellType[rowSeps], sep = ",")

# order x-axis by chrArm, drop: "barcode","cellType","color"
# heatmapMatrix <- heatmapMatrix[,c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21p","21q","22q","Xp","Xq")] # ALL GENES
heatmapMatrix <- heatmapMatrix[,c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q","Xp","Xq")] # HOUSEKEEPING GENES

heatmapMatrix<- as.matrix(x = heatmapMatrix)

# create heatmap
png(filename = paste(myDir,"heatmaps/",patient,".heatmap.NonEpi_epiCycl.housekeeping.png", sep=""), width = 3000, height = 1500)
heatmap.2(x = heatmapMatrix,
          Rowv = FALSE, 
          Colv = FALSE,
          dendrogram = "none", 
          scale = "none",
          trace = "none",
          labRow = groupLabels,
          RowSideColors = as.character(heatmapColors),
          rowsep = rowSeps,
          main = info,
          xlab = "Chromosome Arm", ylab = "Cells",
          col = colorRampPalette(c("red", "white", "blue"))(n = 299), 
          breaks = c(seq(-3.00,-1.01,length=100),
                     seq(-1.00, 1.00,length=100),
                     seq( 1.01, 3.00,length=100)), 
          symbreaks = TRUE,
          symm = FALSE, 
          symkey = TRUE,
          keysize = 0.5,
          cexRow = 1
)
dev.off()
######################


# normalize values on all non-epithelial cells, use hclust to create heatmap clustered on cell values
#######################

info <- paste(patient, "total fraction of UMI counts per cell, centered/scaled by chrArm on all non-epithelial cells, hclust", sep = "; ")

# convert data to wide matrix
heatmapMatrix <- chrArmFracSums[,c("chrArm", "barcode", "UMInorm")] 
heatmapMatrix <- dcast(heatmapMatrix, barcode~chrArm, value.var="UMInorm") 

# fix rownames
rownames(heatmapMatrix) <- heatmapMatrix[,"barcode"] 
# drop barcode columns, reorder by chrArm
# heatmapMatrix <- heatmapMatrix[,c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21p","21q","22q","Xp","Xq")] # ALL GENES
heatmapMatrix <- heatmapMatrix[,c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q","Xp","Xq")] # HOUSEKEEPING GENES

# compute hclust
hclustOrder <- hclust(d = dist(x = heatmapMatrix), method = "average")

# figure out color labels by minor cell type
heatmapColors <- cluster_data.csv$minorCellTypeColor[match(hclustOrder$labels, cluster_data.csv$barcode, nomatch = "NA")]
heatmapColors[heatmapColors==""] <- NA


# add cell type, clust num, and color data
#heatmapMatrix$cellType <- cluster_data.csv$minorCellType[match(heatmapMatrix$barcode, cluster_data.csv$barcode, nomatch = "NA")]
#heatmapMatrix$clustNum <- cluster_data.csv$clustNum[match(heatmapMatrix$barcode, cluster_data.csv$barcode, nomatch = "NA")]
# order y-axis (cells) by cellType, then clustNum
#heatmapMatrix <- heatmapMatrix[order(heatmapMatrix$cellType, heatmapMatrix$clustNum),]
# create row seperators based on cluster numbers
#rowSeps <- c(1,1+which(diff(heatmapMatrix$clustNum)!=0))

# make group row labels
#groupLabels <- rep("", nrow(heatmapMatrix))
#groupLabels[rowSeps] <- paste(heatmapMatrix$clustNum[rowSeps], heatmapMatrix$cellType[rowSeps], sep = ",")

heatmapMatrix<- as.matrix(x = heatmapMatrix)

# create heatmap
png(filename = paste(myDir,"heatmaps/",patient,".heatmap.NonEpi.housekeeping.hclust.png", sep=""), width = 3000, height = 1500)
heatmap.2(x = heatmapMatrix,
          # dendrogram control 
          Rowv = TRUE,
          Colv = FALSE,
          distfun = dist,
          hclustfun = hclust,
          dendrogram = "row", 
          symm = FALSE, 
          
          # data scaling
          scale = "none",
          
          # mapping data to colors
          breaks = c(seq(-3.00,-1.01,length=100),
                     seq(-1.00, 1.00,length=100),
                     seq( 1.01, 3.00,length=100)), 
          symbreaks = TRUE,
          
          # colors
          col = colorRampPalette(c("red", "white", "blue"))(n = 299), 
          
          # block separation 
          rowsep = FALSE, #rowSeps
          
          # level trace
          trace = "none",
          
          # Row/Columns Labeling
          labRow = NULL, #groupLabels
          cexRow = 1,
          # RowSideColors = as.character(heatmapColors),
          
          # color key + density info
          keysize = 0.5,
          symkey = TRUE,
          
          # plot labels
          main = info,
          xlab = "Chromosome Arm", 
          ylab = "Cells"
)
dev.off()
#######################


# old code for finding out how many (housekeeping) genes exist on each chromosome arm
#######################
# old code for finding out how many (housekeeping) genes exist on each chromosome arm
# how many housekeeping genes on each chrArm?
#tmp <- geneData[  geneData$geneName %in% housekeepingGenes  ,  ]
#tmp$chrArm <- paste(tmp$chr, tmp$arm, sep="")
#counts <- aggregate( geneName ~ chrArm  , data = tmp , FUN = function(x){NROW(x)})
#counts
# genesOnChrArm <- read_excel(paste(myDir,"numGenesOnChrs.xlsx",sep = ""))
# ggplot(data=genesOnChrArm, aes(x=chrArm, y=numGenes, group=1)) +
#   geom_line() +
#   geom_point()   
#######################


