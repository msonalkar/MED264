# Purpose: Perform a copy number variation (CNV) analysis on the thymoma data
# Author: Samuel Reynolds
# Date: 11/5/24

library(minfi)
library(conumee)
library(IlluminaHumanMethylationEPICanno.ilm10b3.hg19)
library(dplyr)

# Documentation:
# https://bioconductor.org/packages/release/bioc/vignettes/conumee/inst/doc/conumee.html

# How to combine array types:
# https://support.bioconductor.org/p/97909/

load("~/Desktop/Christensen Lab/Thymoma Classification/Data/combinedPheno.RData")
normal_thymus_samples <- rownames(combinedPheno[combinedPheno$CancerSubtype == "normal thymus",])

controlRownamesMsetGEO <- c("GSM6751733_205555380015_R05C01", "GSM6751734_205555380015_R06C01", "GSM6751735_205555380015_R07C01",
                            "GSM6751736_205555380015_R08C01", "GSM6751737_205555380022_R01C01","GSM6751738_205555380022_R02C01",
                            "GSM6751739_205555380022_R03C01", "GSM6751740_205555380022_R04C01", "GSM6751741_205555380022_R05C01")

#===============================================================================================================
#===============================================================================================================
#  Controls - from GEO
#===============================================================================================================
#===============================================================================================================
# Reading IDATs
setwd("~/Desktop/Christensen Lab/Thymoma Classification/Data/GEO/GSE218549_RAW")
rgSetControls <- read.metharray(basenames = controlRownamesMsetGEO)

# Creating Mset object
MsetControls <- preprocessNoob(rgSetControls)
MsetControls <- convertArray(MsetControls, outType = "IlluminaHumanMethylation450k")

# CNV Intensity Values
controls.intensity.data <- CNV.load(MsetControls)

# Creating CNV annotation object
anno <- CNV.create_anno(bin_minprobes = 10, bin_minsize = 100000, array_type = "overlap", 
                        chrXY = FALSE, exclude_regions = NULL, detail_regions = NULL)

#===============================================================================================================
#===============================================================================================================
#  GEO Thymus Tumor Data
#===============================================================================================================
#===============================================================================================================
# Reading IDATs
setwd("~/Desktop/Christensen Lab/Thymoma Classification/Data/GEO/GSE218549_RAW")
basenamesGEO <- unique(substr(list.files("~/Desktop/Christensen Lab/Thymoma Classification/Data/GEO/GSE218549_RAW"),1,30))
basenamesGEO <- basenamesGEO[!(basenamesGEO %in% controlRownamesMsetGEO)]
rgSetGEO <- read.metharray(basenames = basenamesGEO)

#===============================================================================================================
#===============================================================================================================
#  TGCA Thymus Tumor Data
#===============================================================================================================
#===============================================================================================================
# Reading IDATs
setwd("~/Desktop/Christensen Lab/Thymoma Classification/Data/TCGA/IDATs/TCGA IDATs")
basenamesTCGA <- unique(substr(list.files("~/Desktop/Christensen Lab/Thymoma Classification/Data/TCGA/IDATs/TCGA IDATs"),1,41))
rgSetTCGA <- read.metharray(basenames = basenamesTCGA)

#===============================================================================================================
#===============================================================================================================
#  Combining and normalizing data
#===============================================================================================================
#===============================================================================================================

# USE THIS TO COMBINE
MsetCombined <- combineArrays(object1 = rgSetGEO, object2 = rgSetTCGA, outType = "IlluminaHumanMethylation450k")
MsetCombined <- preprocessNoob(MsetCombined)

#===============================================================================================================
#===============================================================================================================
#  CNV Analysis
#===============================================================================================================
#===============================================================================================================
# CNV Intensity Values
combined.intensity.data <- CNV.load(MsetCombined)

samples <- colnames(combined.intensity.data@intensity)
controls <- colnames(controls.intensity.data@intensity)

##Perform CNV Analysis on ALL SAMPLES
require(lattice)
setwd("~/Desktop/Christensen Lab/Thymoma Classification/CNV Results/GEO Only/Plots")

it = 1
CNVData_List <- list()
for (i in samples){
  print(it)
  it <- it + 1
  x <- CNV.fit(combined.intensity.data[(i)], ref = controls.intensity.data[controls], anno)
  x <- CNV.bin(x)
  x <- CNV.detail(x)
  assign(paste(i,"_CNVData", sep = ""), CNV.segment(x))
  CNVData_List <- append(CNVData_List, get(paste(i,"_CNVData", sep = "")))
  pdf(paste("CNVPlot_", i, ".pdf", sep = ""), height = 9, width = 18)
  print(CNV.genomeplot(get(paste(i,"_CNVData", sep = ""))))
  dev.off()
}

save(CNVData_List, file = "~/Desktop/Christensen Lab/Thymoma Classification/CNV Results/GEO Only/GEO_CNVData_List.RData")

# Creating a dataframe of the bin CNV ratios for all of the samples
combinedCNVRatios <- data.frame()
i = 1
for (dataset in unlist(CNVData_List)){
  print(i)
  i <- i + 1
  name <- dataset@name
  ratios <- as.data.frame(dataset@bin$ratio)
  colnames(ratios) <- name
  ratios <- t(ratios)
  ratios <- as.data.frame(ratios)
  combinedCNVRatios <- rbind(combinedCNVRatios, ratios)
}

combinedPheno_tumorOnly <- read_csv("~/Desktop/Christensen Lab/Thymoma Classification/Data/combinedPheno_tumorOnly.csv")
row.names(combinedCNVRatios)[1:113] <- substr(row.names(combinedCNVRatios)[1:113],1,10)
row.names(combinedCNVRatios)[114:209] <- sub("_noid","",row.names(combinedCNVRatios)[114:209])
combinedCNVRatios <- combinedCNVRatios[row.names(combinedCNVRatios) %in% combinedPheno_tumorOnly$...1,]

View(combinedCNVRatios)
write.csv(combinedCNVRatios, file = "~/Desktop/Christensen Lab/Thymoma Classification/CNV Results/combinedCNVRatios.csv")
save(combinedCNVRatios, file = "~/Desktop/Christensen Lab/Thymoma Classification/CNV Results/combinedCNVRatios.RData")

# Correcting for batch effects

batch <- factor(rep(c("A", "B"), c(106, 94)))
combinedCNVRatios <- t(combinedCNVRatios)
dim(combinedCNVRatios)

all_CNV_batch_corrected <- removeBatchEffect(x=combinedCNVRatios, batch=batch)
all_CNV_batch_corrected <- t(all_CNV_batch_corrected)
write.csv(all_CNV_batch_corrected, file = "~/Desktop/Christensen Lab/Thymoma Classification/CNV Results/all_CNV_batch_corrected.csv")



load("~/Desktop/Christensen Lab/Thymoma Classification/Data/combinedBetas_batch_corrected.RData")
combinedBetas_batch_corrected <- t(combinedBetas_batch_corrected)
combinedBetas_batch_corrected <- combinedBetas_batch_corrected[row.names(combinedBetas_batch_corrected) %in% combinedPheno[,],]
colnames(combinedBetas_batch_corrected) == colnames(all_CNV_batch_corrected)
dim(combinedBetas_batch_corrected)
dim(all_CNV_batch_corrected)









