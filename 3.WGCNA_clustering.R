# WGCNA ATAC-seq only
# ATAC = raw counts filytered 10 CPM

Power = 12
Size = 120
deepsplit = 2
method = "signed hybrid"
remove_outliers = T  # remove outliers? T / F

## load libraries

library(edgeR)
library("DESeq2")
library(RDAVIDWebService)
library(ggplot2)
library("WGCNA")
allowWGCNAThreads()

# directories
input = "../data/"
output = "../results/"

## load data
load(file = paste(input, "dge_atac-seq_10CPM_trim.xz", sep = ""),
     verbose = TRUE)

if (remove_outliers) {
  keep = c(
    "iPSC-sbad2.1",
    "iPSC-sbad3.1",
    "iPSC-neo1.1",
    "DE-sbad2.1",
    "DE-sbad3.1",
    "DE-neo1.1",
    "PGT-sbad2.1",
    "PGT-sbad3.1",
    "PGT-neo1.1",
    "PFG-sbad2.1",
    "PFG-sbad3.1",
    "PFG-neo1.1",
    "PE-sbad2.1",
    "PE-sbad3.1",
    "PE-neo1.1",
    "EP-sbad3.1",
    "EP-neo1.1",
    "EN-sbad3.1",
    "EN-neo1.1",
    "BLC-sbad3.1",
    "BLC-neo1.1"
  )
  dge <-
    dge[, keep, keep.lib.sizes = TRUE]
}

dim(dge)
#174,379     24

## renaming
# ATAC
counts_ATAC = dge$counts


rownames(counts_ATAC)[which(duplicated(rownames(counts_ATAC)))]
# no duplicates


subject = sapply(strsplit(colnames(counts_ATAC), split = "-"), function(x)
  x[2]) # takes second element of column names after splitting by "_"
design = data.frame(
  row.names = colnames(counts_ATAC),
  stage = dge$samples$group,
  subject = subject
) # design is the same for both
dds <- DESeqDataSetFromMatrix(
  countData = counts_ATAC,
  colData = design,
  design = ~ subject + stage
)

design$stage = factor(design$stage, levels(design$stage))

############ Run Vst normalization on the whole dataset (with DESeq2) ############

vsd <- varianceStabilizingTransformation(dds)
vstMat = assay(vsd)
write.table(
  vstMat,
  file = paste(output, "vstMat.txt", sep = ""),
  sep = "\t",
  quote = F
)

#### blind=TRUE should be used for comparing samples in an manner unbiased by prior  ###
#### information on samples, for example to perform sample QA (quality assurance).   ###
#### blind=FALSE should be used for transforming data for downstream analysis, where ###
#### the full use of the design information should be made							 ###

vsd_b <- varianceStabilizingTransformation(dds, blind = F)
vstMat_b = assay(vsd)
write.table(
  vstMat_b,
  file = paste(output, "vstMat.blind_F.txt", sep = ""),
  sep = "\t",
  quote = F
)

### MDS plot
pdf(paste(output, "mds.pdf", sep = ""), width = 10)
par(mfrow = c(1, 2))
plotMDS(vstMat, col = as.numeric(dds$stage))
plotMDS(vstMat_b, col = as.numeric(dds$stage))
dev.off()
## no difference in vst transform here


##################### run WGCNA  #####################
if (!file.exists("net.xz")) {
  degData = t(vstMat)
  gsg = goodSamplesGenes(degData, verbose = 3)
  
  gsg$allOK
  
  ### pick power
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft = pickSoftThreshold(degData, powerVector = powers, verbose = 5)
  pdf(
    file = paste(output, "WGCNA_softThreshold.pdf", sep = ""),
    width = 9,
    height = 5
  )
  
  par(mfrow = c(1, 2))
  
  cex1 = 0.9
  
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence")
  )
  
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.90, col = "red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  dev.off()
  
  # looks like power of 12 is fine
  ##### run WGCNA module detection
  
  net = blockwiseModules(
    degData,
    power = Power,
    maxBlockSize = 15500,
    TOMType = "signed",
    networkType = method,
    minModuleSize = Size,
    mergeCutHeight = 0.10,
    detectCutHeight = 0.99,
    deepSplit = deepsplit,
    numericLabels = TRUE,
    pamStage = T,
    pamRespectsDendro = T,
    loadTOM = T,
    verbose = 3
  )
  save(net,
       file = paste(output, "net.xz", sep = ""),
       compress = "xz")
  table(net$colors)
  
  
  pdf(
    paste(output, "WGCNA_dendrogram.30M.pdf", sep = ""),
    width = 12,
    height = 9
  )
  moduleColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    net$dendrograms[[1]],
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  
  
  
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs
  
  geneTree = net$dendrograms[[1]]
  
  
  save(
    MEs,
    moduleLabels,
    moduleColors,
    geneTree,
    file = paste(output, "ATAC_only_network.RData", sep = "")
  )
  
  lsdatME = moduleEigengenes(degData, moduleColors)$eigengenes
  datME = moduleEigengenes(degData, moduleColors)$eigengenes
  dissimME = (1 - t(WGCNA::cor(datME, method = "p"))) / 2
  hclustdatME = hclust(as.dist(dissimME), method = "average")
  
  MEs = moduleEigengenes(degData, moduleColors)$eigengenes
  pheno = data.frame(stage = as.numeric(design$stage))
  
  MET = orderMEs(cbind(MEs, pheno))
  pdf(
    paste(output, "WGCNA_cluster_eigenvectors.30M.pdf", sep = ""),
    width = 8,
    height = 12
  )
  par(cex = 0.9)
  plotEigengeneNetworks(
    MET,
    "",
    marDendro = c(0, 4, 1, 2),
    marHeatmap = c(3, 4, 1, 2),
    cex.lab = 0.8,
    xLabelsAngle = 90
  )
  dev.off()
  
  
  pdf(
    paste(output, "Module_eigengenes.30M.barplots.pdf", sep = ""),
    width = 10,
    height = 3
  )
  for (i in 1:length(unique(moduleColors))) {
    which.module = unique(moduleColors)[i]
    ME = datME[, paste("ME", which.module, sep = "")]
    barplot(
      ME,
      col = as.numeric(dds$stage),
      ylab = "eigengene expression",
      xlab = "array sample",
      main = which.module
    )
  }
  dev.off()
  
  gene2module = data.frame(gene = colnames(degData), module = moduleColors)
  write.table(
    gene2module,
    file = paste(output, "WGCNA.gene2module.30M.txt", sep = ""),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  #### make ribbon plots for each module  #####
  pdf(paste(output, "Modules.ribbon_plots.pdf", sep = ""))
  for (j in 1:dim(datME)[2]) {
    stages = levels(design$stage)
    min = rep(0, length(stages))
    max = rep(0, length(stages))
    mean = rep(0, length(stages))
    for (i in 1:length(stages)) {
      min[i] = min(datME[which(design$stage == stages[i]), j])
      max[i] = max(datME[which(design$stage == stages[i]), j])
      mean[i] = mean(datME[which(design$stage == stages[i]), j])
    }
    plot_data = data.frame(
      stages = c(1:length(stages)),
      min = min,
      max = max,
      mean = mean
    )
    
    p = ggplot(plot_data, aes(stages, mean)) + geom_ribbon(aes(ymin = min, ymax =
                                                                 max),
                                                           colour = "lightgrey",
                                                           fill = "lightgrey") + geom_line(color = "steelblue4", lwd =
                                                                                             1) + theme_bw() + scale_x_continuous(breaks = c(1:8)) + ylab("module eigengene") +
      ggtitle(colnames(datME)[j])
    
    p = ggplot(plot_data, aes(stages, mean)) + geom_ribbon(aes(ymin = min, ymax =
                                                                 max),
                                                           colour = "lightgrey",
                                                           fill = "lightgrey") + geom_line(color = "steelblue4", lwd =
                                                                                             1) + theme_bw() + scale_x_continuous(breaks = c(1:8), labels = stages) +
      ylab("module eigengene") + ggtitle(colnames(datME)[j])
    
    print(p)
  }
  dev.off()
}

if (file.exists("net.xz")) {
  vstMat = read.table(paste(output, "vstMat.txt", sep = ""))
  degData = t(vstMat)
  load(paste(output, "net.xz", sep = ""))
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  ### get connectivity values for each gene within its module
  ### get connectivity values:
  ADJ1 = abs(WGCNA::cor(degData, use = "p")) ^ Power
  Alldegrees1 = intramodularConnectivity(ADJ1, moduleColors)
  Alldegrees1$Module = moduleColors
  
  write.table(
    Alldegrees1,
    file = paste(output, "WGCNA.connectivity.gene2module.txt", sep = ""),
    sep = "\t",
    quote = F
  )
  
  # split table based on modules
  
  Alldegrees1_list <- split(Alldegrees1 , f = Alldegrees1$Module)
  for (i in 1:length(Alldegrees1_list)) {
    write.table(
      Alldegrees1_list[i],
      file = paste0(
        output,
        "WGCNA.connectivity.gene2module.",
        names(Alldegrees1_list)[i],
        ".txt"
      ),
      sep = "\t",
      quote = F
    )
  }
  
}