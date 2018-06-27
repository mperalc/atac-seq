# correlated peaks 1MB either side of peaks overlapping credible sets

library(GenomicRanges)
library(WGCNA)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(gridGraphics)
library(dplyr)
source("/Users/Marta/Documents/WTCHG/R scripts/atac-seq/corWGCNAmodules.R")

options(stringsAsFactors = FALSE)

# 1. correlate ATAC and RNA WGCNA modules

ATAC_dir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/10CPM_P12_S120_deepSplit2_signedHybrid_noOutliers/"
RNA_dir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/WGCNA_P12/"

# Read in module files for RNA and ATAC-seq
ATAC = read.table(file = paste(ATAC_dir, "WGCNA.gene2module.30M.txt", sep =
                                 ""),
                  header = T)
RNA = read.table(file = paste(RNA_dir, "WGCNA.gene2module.30M.txt", sep =
                                ""),
                 header = T)


# read in expression/atac data
# RNA
load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz")  # 15221 genes and lincRNA
# ATAC
load(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq_10CPM_trim.xz"
)

# plotting correlated modules
setwd(ATAC_dir)
corWGCNAmodules(modulesDataSet1 = ATAC, modulesDataSet2 = RNA, dataSet1 = dge, dataSet2 = dge_cc)

# normalize expression data
datExpr = cbind(dge_cc$genes[1:2], dge_cc$counts)
datExpr[3:ncol(datExpr)] = cpm(datExpr[3:ncol(datExpr)])

# Normalize ATAC data
datATAC = cbind(dge$genes[1:2], dge$counts)
datATAC[3:ncol(datATAC)]  = cpm(datATAC[3:ncol(datATAC)])

############# look 500kb around every credible set ##########
### check loci classification from Anubha - there are some SNPs shared between loci

credset_names = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/credset_names.txt",
                           header = F)

credset = list()

for (n in credset_names$V1) {
  credset[[n]] = read.table(
    file = paste(
      "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credset/",
      n,
      sep = ""
    ),
    header = T
  ) # 380 99% credible sets T2D 2017 HRC (Mahajan et al 2018)
  
}

# peaks over credible sets
peaksOverSNPs = read.table(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/over10CPM_peaks_over_SNPs.txt"
)
nrow(peaksOverSNPs)

# allocate each peak in correct credible set
# some will be empty, others will have peak
# integrate peak info and credset SNPs with findOverlaps


peaksOverSNPs.gr  = makeGRangesFromDataFrame(peaksOverSNPs)

peaks_per_credset = list()
for (n in names(credset)) {
  credset_subset = credset[[n]]
  credset_subset = credset_subset[, c("Chr", "Pos")]
  colnames(credset_subset) = c("chr", "start")
  credset_subset$end = credset_subset$start + 1
  credset_subset$chr = paste("chr", credset_subset$chr, sep = "")
  credset_subset$chr = as.character(credset_subset$chr)
  # duplicated within same credible set loci??
  credset_subset = unique(credset_subset)
  rownames(credset_subset) = paste("snp",
                                   credset_subset$chr,
                                   credset_subset$start,
                                   credset_subset$end,
                                   sep = "_")
  credset.gr = makeGRangesFromDataFrame(credset_subset)
  
  # Get all peaks in module overlapping credible set SNPs :
  overlapPeak = subsetByOverlaps(peaksOverSNPs.gr , credset.gr)
  hits <-
    findOverlaps(peaksOverSNPs.gr, credset.gr) # to add SNP metadata to overlapPeak table
  rsid <-
    CharacterList(split(names(credset.gr)[subjectHits(hits)], queryHits(hits)))
  rsid = paste(rsid, collapse = ",")
  
  mcols(overlapPeak) <- DataFrame(mcols(overlapPeak), rsid)
  
  
  overlapPeak = as.data.frame(overlapPeak)
  overlapPeak = overlapPeak[, c(1:3, 6)]
  colnames(overlapPeak) = c("Chr", "Start", "End", "SNPpos")
  if (nrow(overlapPeak) > 0) {
    overlapPeak$PeakID = paste("peak", overlapPeak$Chr, overlapPeak$Start, sep = "_")
  }
  peaks_per_credset[[n]] = overlapPeak
  
}
rm(credset_subset, overlapPeak)

# how many credsets don't have any overlapping peaks?
lapply(peaks_per_credset, nrow)

# number of credible sets w/o any peaks over their SNPs
peaks_per_credset %>%
  lapply(nrow) %>%
  (function(x)
    x == 0) %>%
  sum

nrow(credset_names) # number credible sets
# 117 of 380 credible sets don't have any overlapping peaks (check this number every time)
(177 / 380) * 100 # 47%

# 53% do have
peaks_per_credset %>%
  lapply(nrow) %>%
  unlist %>%
  sum
# 1229 peaks in total in all credible sets.
nrow(peaksOverSNPs) # 1007? 222 duplicated?
temp = do.call("rbind", peaks_per_credset)
length(unique(temp$PeakID)) # 1007. So there are duplicates
sum(duplicated(temp$PeakID))
dups = temp[duplicated(temp$PeakID), "PeakID"] # list of duplicated ids

temp = temp[which(temp$PeakID %in% dups),]


#### get PPA for all SNPs under peaks in every credible set
# list of credible set
# SNP overlap peak
# snp, ppa, peak
# if no peak overlapping - "none" under peak col
# cbind and save with row names
SNPs_in_peaks = list()

for (n in names(credset)) {
  peaks = peaks_per_credset[[n]]
  
  credset_subset = credset[[n]]
  credset_subset = credset_subset[, c("Chr", "Pos")]
  colnames(credset_subset) = c("chr", "start")
  credset_subset$end = credset_subset$start + 1
  credset_subset$chr = paste("chr", credset_subset$chr, sep = "")
  credset_subset$chr = as.character(credset_subset$chr)
  credset_subset = unique(credset_subset)
  rownames(credset_subset) = paste("snp",
                                   credset_subset$chr,
                                   credset_subset$start,
                                   credset_subset$end,
                                   sep = "_")
  credset.gr = makeGRangesFromDataFrame(credset_subset)
  
  # Get credible set SNPs overlapping peaks in module :
  overlapSNP = subsetByOverlaps(credset.gr, peaksOverSNPs.gr)
  hits <-
    findOverlaps(credset.gr, peaksOverSNPs.gr) # to add SNP metadata to overlapPeak table
  peakid <-
    unlist(split(names(peaksOverSNPs.gr)[subjectHits(hits)], queryHits(hits)))
  #peakid = paste(peakid, collapse = ",")
  
  mcols(overlapSNP) <- DataFrame(mcols(overlapSNP), peakid)
  
  overlapSNP = as.data.frame(overlapSNP)
  if (nrow(overlapSNP) == 0) {
    credset[[n]]$IndexSNP = paste(credset[[n]]$Chr, credset[[n]]$Pos, sep =
                                   "_")
    overlapSNP = credset[[n]]
    overlapSNP$PeakID = rep("none", nrow(overlapSNP))
    SNPs_in_peaks[[n]] = overlapSNP
  }
  else {
  overlapSNP = overlapSNP[, c(1:2, 6)]
  colnames(overlapSNP) = c("Chr", "Pos", "PeakID")
  # if (nrow(overlapSNP) > 0) {
  #   overlapSNP$PeakID = paste("peak", overlapSNP$Chr, overlapSNP$Start, sep = "_")
  # }
  
  overlapSNP$IndexSNP = paste(gsub("\\chr*", "", overlapSNP$Chr), overlapSNP$Pos, sep =
                                "_")
  credset[[n]]$IndexSNP = paste(credset[[n]]$Chr, credset[[n]]$Pos, sep =
                                  "_")
  # SNPs_in_peaks[[n]] = credset[[n]]
  # matches = which(credset[[n]]$IndexSNP %in% overlapSNP$IndexSNP)
  # SNPs_in_peaks[[n]]$PeakID = rep("none", nrow(SNPs_in_peaks[[n]]))
  # SNPs_in_peaks[[n]][matches,"PeakID"] <- overlapSNP$PeakID
  SNPs_in_peaks[[n]] = merge(credset[[n]], overlapSNP[,c(3,4)], by = "IndexSNP", all = T)
  SNPs_in_peaks[[n]][is.na(SNPs_in_peaks[[n]]$PeakID),"PeakID"] = "none"
  }
}
SNPs_in_peaks = do.call("rbind",SNPs_in_peaks)
write.table(SNPs_in_peaks,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/credset_SNPs_in_peaks_PPA.txt", 
            sep="\t", quote = F, row.names = T, col.names = T)

######### then get df for each credible set + 500kb each side

credset_plus500 = credset
gene_matrices_dist = list()
gene_matrices_pval = list()
gene_matrices_qval = list()
gene_matrices_cor = list()

library(biomaRt)

for (n in names(credset_plus500)) {
  credset_plus500[[n]] = credset_plus500[[n]][, c("Chr", "Pos")]
  colnames(credset_plus500[[n]]) = c("chr", "start")
  credset_plus500[[n]]$end = credset_plus500[[n]]$start + 1
  chr = unique(credset_plus500[[n]]$chr)
  start = min(credset_plus500[[n]]$start) - 500000
  end =  max(credset_plus500[[n]]$end) + 500000
  
  # filterlist = list(paste(gsub("\\chr*", "", chr) , #remove "chr" in every element in location
  #                         start, end, sep = ":"))
  
  filterlist = list(paste(chr, start, end, sep = ":"))
  
  ensembl = useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    host = "grch37.ensembl.org",
    path = "/biomart/martservice" ,
    dataset = "hsapiens_gene_ensembl"
  )
  # grch37 used for annotation of RNA-seq background
  
  # search distance to genes and correlation with peak for those loci that have overlapping peaks:
  if (nrow(peaks_per_credset[[n]]) > 0) {
    results <-
      getBM(
        attributes = c(
          'ensembl_gene_id',
          'external_gene_name',
          'gene_biotype',
          'chromosome_name',
          "start_position",
          "end_position"
        ),
        filters = c("chromosomal_region"),
        values = filterlist,
        mart = ensembl
      )
    
    
    results = results[which(results$gene_biotype == "lincRNA" |
                              results$gene_biotype == "protein_coding"),]
    
    results = results[!duplicated(results$external_gene_name), ]
    rownames(results) = results$external_gene_name
    results2 = results[, c("chromosome_name", "start_position", "end_position")]
    colnames(results2) = c("chr", "start", "end")
    results2$chr = paste("chr", results2$chr, sep = "")
    results.gr = makeGRangesFromDataFrame(results2)
    mcols(results.gr) <-
      DataFrame(mcols(results.gr),
                results$ensembl_gene_id,
                results$gene_biotype)
    
    # get atac peaks in credset
    peaks.gr = makeGRangesFromDataFrame(peaks_per_credset[[n]][, c(1:3)])
    
    # initialize df
    gene_matrices_dist[[n]] = data.frame(matrix(
      ncol = length(results.gr@ranges@NAMES),
      nrow = peaks.gr@seqnames@lengths
    ))
    colnames(gene_matrices_dist[[n]]) = results.gr@ranges@NAMES
    rownames(gene_matrices_dist[[n]]) = peaks.gr@ranges@NAMES
    
    gene_matrices_pval[[n]] = gene_matrices_dist[[n]] # save for later
    gene_matrices_cor[[n]] = gene_matrices_dist[[n]] # save for later
    # populate df calculating distances between every peak and every gene
    for (p in 1:peaks.gr@seqnames@lengths) {
      gene_matrices_dist[[n]][p, ] = distance(peaks.gr[p], results.gr)
    }
    
    ######################################## calculating correlations
    # subset df for genes present
    #Expr = t(datExpr)
    #datExpr_subset = Expr[which(datExpr$external_gene_name %in% colnames(gene_matrices_dist[[n]])),]
    datExpr_subset = datExpr[which(datExpr$external_gene_name %in% colnames(gene_matrices_dist[[n]])),]
    
    # if there is only one peak:
    if (is.null(nrow(datExpr_subset))) {
      datExpr_subset = as.data.frame(datExpr_subset)
      colnames(datExpr_subset) = dge_cc$genes[dge_cc$genes$external_gene_name %in% colnames(gene_matrices_dist[[n]]), "external_gene_name"]
      datExpr_subset = t(datExpr_subset)
      
    }
    else{
      rownames(datExpr_subset) = dge_cc$genes[dge_cc$genes$ensembl_gene_id %in% rownames(datExpr_subset), "external_gene_name"]
    }
    # subset atac count table for peaks present
   # ATAC = t(datATAC)
   # if (is.null(ncol(ATAC[which(rownames(ATAC) %in% rownames(gene_matrices_dist[[n]])),]))) {
      if (is.null(ncol(datATAC[which(rownames(datATAC) %in% rownames(gene_matrices_dist[[n]])),]))) {
        
      # datATAC_subset = ATAC[which(rownames(ATAC) %in% rownames(gene_matrices_dist[[n]])),]
      # datATAC_subset = as.data.frame(datATAC_subset)
      # cormatrix = as.data.frame(cor(t(datExpr_subset), (datATAC_subset)), use = "p")
      # colnames(cormatrix) = rownames(gene_matrices_dist[[n]])
      # # correcting weird behavior when there is only one peak
      # corAndPvalue(x = t(datExpr_subset), y = datATAC_subset )
        
        datATAC_subset = datATAC[which(rownames(datATAC) %in% rownames(gene_matrices_dist[[n]])),]
        datATAC_subset = as.data.frame(datATAC_subset)
        #cormatrix = as.data.frame(cor(t(datExpr_subset), t(datATAC_subset), use = "p"))
        cormatrix = corAndPvalue(x = t(datExpr_subset[3:ncol(datExpr_subset)]), y = t(datATAC_subset[3:ncol(datATAC_subset)]) )
        pvals = cormatrix$p
        qvals = apply(pvals,2, function(x){p.adjust(x,method = "bonferroni")})
        cormatrix = cormatrix$cor
        colnames(cormatrix) = rownames(gene_matrices_dist[[n]])
        # correcting weird behavior when there is only one peak
      
    }
    else{
      datATAC_subset = datATAC[which(rownames(datATAC) %in% rownames(gene_matrices_dist[[n]])),]
      datATAC_subset = as.data.frame(datATAC_subset)
      
      # calculate pearson correlation between peaks and genes
      #cormatrix = as.data.frame(cor(t(datExpr_subset[3:ncol(datExpr_subset)]),t(datATAC_subset[3:ncol(datATAC_subset)]), use = "p"))
      cormatrix = corAndPvalue(x = t(datExpr_subset[3:ncol(datExpr_subset)]), y = t(datATAC_subset[3:ncol(datATAC_subset)]) )
      pvals = cormatrix$p
      qvals = apply(pvals,2, function(x){p.adjust(x,method = "bonferroni")})
      if (nrow(pvals) == 1) {
        qvals = t(as.data.frame(qvals))
        
        rownames(qvals) = rownames(pvals)
      }
      cormatrix = cormatrix$cor
      
    }
    cormatrix = t(cormatrix)
    pvals = t(pvals)
    qvals = t(qvals)
    
   # pvals = as.data.frame(corPvalueStudent(as.matrix(cormatrix), nSamples))
    # for (c in colnames(cormatrix)) {
    #   gene_matrices_cor[[n]][, c] = cormatrix[, c]
    #   gene_matrices_pval[[n]][, c] = pvals[, c]
    #   #gene_matrices_qval[[n]][, c] = qvals[, c]
    #   
    # }
    gene_matrices_cor[[n]] = cormatrix
    gene_matrices_pval[[n]] = pvals
    gene_matrices_qval[[n]] = qvals
    
    
  }
}


# save matrices
save(
  list = c(
    "gene_matrices_dist",
    "gene_matrices_cor",
    "gene_matrices_pval",
    "gene_matrices_qval"
  ),
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/gene_matrices_dist_cor_pval_qval.RData"
)

load(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/gene_matrices_dist_cor_pval_qval.RData")

# get all genes in that position
# rows: every peak present in that credset
# cols are genes: get p-val of GS with ATAC module
# or just direct pearson correlation? check how both matrices would look
# say NA for genes that are not in dge (<1CPM etc)
# anothr matrix for distance to peak

# get some data: number of sign cor genes, avg distance per peak
# keep in mind significantly correlated genes can be positively or negatively correlated

# number of sign cor genes. per peak Ignore genes not expressed
# signCorGenesPerPeak = lapply(gene_matrices_pval, function(x) {
#   rowSums(x < 0.05, na.rm = T)
# })
# 
# png(filename = paste(ATAC_dir, "hist_gene_cor_per_peak.png", sep = ""))
# hist(
#   unlist(signCorGenesPerPeak),
#   breaks = c(70),
#   main = "How many genes cor. per peak?",
#   xlab = "N genes sign. correlated per peak"
# )
# dev.off()

# for genes within 500kb of every peak
signCorGenesPerPeak_500 = mapply(function(X, Y) {
    rowSums(Y < 500000 & X < 0.05, na.rm = T)
}, X = gene_matrices_pval, Y = gene_matrices_dist)

png(filename = paste(ATAC_dir, "hist_gene_cor_per_peak_500kbfromPeak.png", sep = ""))
hist(
  unlist(signCorGenesPerPeak_500),
  breaks = c(70),
  main = "How many genes cor. per peak?",
  xlab = "N genes sign. correlated per peak"
)
dev.off()

table(unlist(signCorGenesPerPeak_500))
summary(unlist(signCorGenesPerPeak_500))
# per loci
# lapply(gene_matrices_pval, function(x) {
#   hist(rowSums(x < 0.05, na.rm = T))
# })

# how many genes sign. positively cor per peak?

AND1 <- function(...)
  Reduce("&", list(...))

combined = list()
for (n in names(gene_matrices_pval)) {
  for (r in 1:nrow(gene_matrices_pval[[n]])) {
    a1 = gene_matrices_pval[[n]][r, ] < 0.05
    a2 = gene_matrices_cor[[n]][r, ] > 0
    combined[[n]] = rbind(combined[[n]], AND1(a1, a2))
  }
}

signPosCorGenesPerPeak = lapply(combined, function(x) {
  rowSums(x, na.rm = T)
})


png(filename = paste(ATAC_dir, "hist_gene_pos_cor_per_peak.png", sep = ""))
hist(
  unlist(signPosCorGenesPerPeak),
  breaks = c(70),
  main = "How many genes cor. (r>0) per peak?",
  xlab = "N genes sign. correlated (r>0) per peak"
)
dev.off()
table(unlist(signPosCorGenesPerPeak))
summary(unlist(signPosCorGenesPerPeak))

#  distance to closest genes per peak

distToGenes = lapply(gene_matrices_dist, function(x) {
  apply(x, 1, function(y)
    min(y))
  
})
png(filename = paste(ATAC_dir, "dist_to_closest_gene_per_peak.png", sep = ""))
hist(
  unlist(distToGenes),
  breaks = c(10000),
  xlim = c(0, 1000),
  main = "How far is each peak from closest gene?",
  xlab = "Distance from closest gene (bp)"
)
dev.off()
table(unlist(distToGenes))
summary(unlist(distToGenes))


# distance to closest significantly correlated peaks


combined = list()
for (n in names(gene_matrices_pval)) {
  v = list()
  for (r in 1:nrow(gene_matrices_pval[[n]])) {
    a1 = gene_matrices_pval[[n]][r, ] < 0.05
    a2 = gene_matrices_cor[[n]][r, ] > 0
    temp = gene_matrices_dist[[n]][r, which(AND1(a1, a2))]
    if (is.integer(temp)) {
      names(temp) = rownames(gene_matrices_dist[[n]][r, ])
      v[[r]] = temp
    }
    
  }
  combined[[n]] = v
}

png(filename = paste(ATAC_dir, "dist_to_closest_sign_cor_gene_per_peak.png", sep = ""))
hist(
  unlist(combined),
  breaks = c(100),
  main = "How far is each peak from closest sign. pos. cor. gene?",
  xlab = "Distance from closest sign. pos. cor. gene (bp)"
)
dev.off()
table(unlist(combined))
summary(unlist(combined))
sum(table(unlist(combined)))


# distance to most significantly correlated peak


# check peaks with most correlations

# EYA2 - SLC2A10 good candidate
gene_matrices_pval$credible_set_Eur_EYA2_20_45317678.txt["peak_chr20_45317502",]
gene_matrices_cor$credible_set_Eur_EYA2_20_45317678.txt["peak_chr20_45317502",]
gene_matrices_dist$credible_set_Eur_EYA2_20_45317678.txt["peak_chr20_45317502",]

# KCNQ1 - PPA 0.81263 nothing cor and also expressed

gene_matrices_pval$credible_set_Eur_KCNQ1_11_2858546.txt["peak_chr11_2858278",]
gene_matrices_cor$credible_set_Eur_KCNQ1_11_2858546.txt["peak_chr11_2858278",]
gene_matrices_dist$credible_set_Eur_KCNQ1_11_2858546.txt["peak_chr11_2858278",]

# PROX1 - PPA 0.64 nothing pos. correlated
gene_matrices_pval$credible_set_Eur_PROX1_1_214150821.txt["peak_chr1_214150098",]
gene_matrices_cor$credible_set_Eur_PROX1_1_214150821.txt["peak_chr1_214150098",]
gene_matrices_dist$credible_set_Eur_PROX1_1_214150821.txt["peak_chr1_214150098",]

# SLC30A8- PPA 0.35 nothing pos. correlated

gene_matrices_pval$credible_set_Eur_SLC30A8_8_118185025.txt["peak_chr8_118183924",]
gene_matrices_cor$credible_set_Eur_SLC30A8_8_118185025.txt["peak_chr8_118183924",]
gene_matrices_dist$credible_set_Eur_SLC30A8_8_118185025.txt["peak_chr8_118183924",]

# ATP1B2 - PPA 0.34 nothing interesting

gene_matrices_pval$credible_set_Eur_ATP1B2_17_7740170.txt["peak_chr17_7739692",c(29:95)]
gene_matrices_cor$credible_set_Eur_ATP1B2_17_7740170.txt["peak_chr17_7739692",c(29:95)]
gene_matrices_dist$credible_set_Eur_ATP1B2_17_7740170.txt["peak_chr17_7739692",c(29:95)]


# IGF2BP3- PPA 0.33 strong negative cor with IGF2BP3
# pos cor with CCDC126, 0.60 r with 0.001926281 pval

gene_matrices_pval$credible_set_Eur_IGF2BP3_7_23434606.txt["peak_chr7_23433911",]
gene_matrices_cor$credible_set_Eur_IGF2BP3_7_23434606.txt["peak_chr7_23433911",]
gene_matrices_dist$credible_set_Eur_IGF2BP3_7_23434606.txt["peak_chr7_23433911",]






# are significant peaks significantly closer than the overall population of peaks?
# sum of significant cor per peak: any peak sums most significant correlations??
# eg. peak_chr17_3862531

