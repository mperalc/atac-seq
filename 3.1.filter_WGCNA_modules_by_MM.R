# Filter elements within WGCNA modules by module membership
# rename 3.2
# Relating ATAC peak modules of WGCNA to the modules of RNA-seq

library(WGCNA)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(gridGraphics)
library(modeest)

options(stringsAsFactors = FALSE)

# locations
ATAC_dir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/WGCNA_P12_S120_deepSplit2_signedHybrid_rmOutliersT/"

# Read in module files for RNA and ATAC-seq
ATAC = read.table(file = paste(ATAC_dir, "WGCNA.gene2module.30M.txt", sep = ""),
                  header = T)

# read in expression/atac data
# ATAC
load(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq_10CPM_trim_conservativeCounts_allQualitySamples.xz"
)

# Number of ATAC-seq peaks, genes and samples
nPeaks = nrow(ATAC)
nSamples = 24

datATAC = t(dge$counts) # ATAC data needs to be in format row-sample and col-peak
# The eigengene is the summary profile of each module we are going to correlate

MEs0 = moduleEigengenes(datATAC, ATAC$module)$eigengenes
MEs = orderMEs(MEs0)


# correlation of eigenvectors - results of WGCNA
moduleTraitCor = cor(MEs, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# df to examine more easily
corDF = as.data.frame(moduleTraitCor)

########### graphical representation of the correlation
sizeGrWindow(10, 6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2),
                   "\n(",
                   signif(moduleTraitPvalue, 1),
                   ")",
                   sep = "")
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))

tiff(
  file = paste(ATAC_dir, "Heatmap_cor_ATAC_modules_WGCNA.tiff", sep = ""),
  width = 25,
  height = 25,
  units = "in",
  compression = "lzw",
  res = 300
)
# Display the correlation values within a heatmap plot
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(MEs),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = greenWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)

dev.off()


# For each module, we define a quantitative measure of module membership MM as the correlation of the module eigengene 
# and the peak expression profile. This allows us to quantify the similarity of all peaks on the array to every module

# names (colors) of the modules for ATAC
modNames = substring(names(MEs), 3)

# calculate MM
peakModuleMembership = as.data.frame(cor(datATAC, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(peakModuleMembership), nSamples))
names(peakModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

# get for every module a scatterplot with cor (x) and p-value(y) for all peaks 
# from this dynamically stablish where to put the line for inclusion in module?

plot_list1 = list()
to_inspect = list()

for (m in names(peakModuleMembership)) {
      module = substring(m, 3)
      moduleATAC = ATAC$module == module
      column = match(m,colnames(peakModuleMembership))
      sizeGrWindow(7, 7)
      par(mfrow = c(1, 1))
      
      # # Save non-ggplot plot to an object using a null PDF device
      # # http://stackoverflow.com/a/14742001/120898
      pdf(NULL)
      dev.control(displaylist = "enable")
      
      # create data frame to fit
       x = abs(peakModuleMembership[moduleATAC, column])
       y = abs(MMPvalue[moduleATAC,column])
       
       # to inspect afterward
       to_inspect[[module]] = as.data.frame(cbind(rownames(MMPvalue)[moduleATAC],MMPvalue[moduleATAC,column]))
       colnames(to_inspect[[module]]) = c("peak","p")
       to_inspect[[module]][["p"]] = as.numeric(unlist(to_inspect[[module]]["p"]))
       
       mode = mlv(y, method = "mfv")
       to_inspect[[module]]$mode = rep(mode$M,nrow(to_inspect[[module]]))
       
      # df = data.frame(x,y)
      # model = nls(y ~ SSasymp(x,Asym,r0,lrc), data = df, trace = T, control = nls.control(maxiter = 500)) 
      # can't fit nls with data that has zero-residual !! (my data)
       # workaround with mode of the p-value, that marks a value similar to the asymptote
       # peaks to ditch are coloured in red
      plot(x, y,
           col = ifelse(y > mode$M, "red", "black"),
           xlab = paste("Module Membership in", module, "module"),
        ylab = paste(
          "MM p-value for",
          substring(m, 3),
          "module",
          sep = " "
        ),
        main = paste("Module membership cor vs. p-value in", module, "module \n")
      )
     
      plot_list1[[m]] <- recordPlot()
      invisible(dev.off())
      
      plot_list1[[m]] = plot_to_gtable(plot_list1[[m]])  # save to gTable to arrange them later
    }
  
  # plot list from above
  # create progress bar
library(tcltk)
total = length(names(peakModuleMembership))
pb <- tkProgressBar(
  title = "progress bar",
  min = 0,
  max = total,
  width = 300
)

output_name = "scatterplot_MM"
# I have to open pdf connection before looping so that it saves one gene on each page
somePDFPath = paste(ATAC_dir,
                    output_name,
                    ".pdf",
                    sep = "")
pdf(file = somePDFPath)
j = 1
for (i in names(peakModuleMembership)) {
  if (identical(plot_list1[[i]], NULL)) {
    next
  }
  else{
    plot(plot_list1[[i]])
  }
  Sys.sleep(0.1)
  setTkProgressBar(pb, j, label = paste(round(j / total * 100, 0),
                                        "% done"))
  j = j + 1
}
close(pb)
dev.off()

######### inspect some elements removed from the modules above

to_inspect = lapply(to_inspect, function(df){
  df[order(df$p,decreasing = F),]
})

#### save WGCNA file with peaks over X p-value
ATAC_reduced = do.call(rbind, to_inspect)
ATAC_reduced = ATAC_reduced[ATAC_reduced$p < ATAC_reduced$mode, ]
ATAC = ATAC[ATAC$gene %in% ATAC_reduced$peak, ]

write.table(ATAC, file = paste(ATAC_dir, "ATAC_pval_selected_MM.txt", sep = ""), 
            sep = "\t", quote = F, row.names = F)  
save(to_inspect, file = paste(ATAC_dir, "ATAC_list_all_MMpvals.Rdata", sep = ""))  # these are all peaks with pval of MM for each module

