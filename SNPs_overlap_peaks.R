# overlap ATAC peaks and T2D GWAS SNPs

### load libraries
library(dplyr)
#library(biomaRt)
#library(readr)
#library(vcfR)
library(GenomicRanges)
library(BSgenome) # for SNP by overlaps function

#### input files
counts <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/counts_all_samples_peaks_narrowpeaks_conservative_set_highQualityBams_renamed.txt",
                     header = T) # bed file with start and end of peaks
peaks = counts[1:3]
rm(counts)

credset = read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/DIAGRAM_T2D_1000G_99_Credsets/credible_sets_diagram_forMAGENTA.txt",
                     header = F) # 99 credible sets T2D 2017 DIAGRAM
credset=unique(credset[1:2]) # SNPs with unique positions (some might be shared between loci?)
colnames(credset)=c("Chr","Pos")
nrow(peaks)
any(duplicated(peaks)) # there are no duplicated peaks
#[1] FALSE
nrow(credset)

# save the step of assigning SNP id for after intersecting peaks
# start_time <- Sys.time()
# snps = read.vcfR( "/Users/Marta/Documents/WTCHG/DPhil/Data/SNPs_b150_GRCh37p13/common_all_20170710.vcf.gz", verbose = FALSE ) # common SNPs build 150 GRCh37p13
# end_time - start_time

# remove sex chromosomes
r1=rownames(peaks[grep("chrY", rownames(peaks), value=F),])
r2=rownames(peaks[grep("chrX", rownames(peaks), value=F),])

peaks=peaks[!rownames(peaks) %in% c(r1,r2),]
rm(r1,r2)
nrow(peaks)

########### 
# 
# # check SNP names
# # load for SNPs 
# ensembl_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_snp")
# 
# 
# getSNPsInDf = function(credset){
#   # biomart query forward strand
#   filterlist <- list(paste(credset[,"Chr"],credset[,"Pos"],credset[,"Pos"],1,sep = ":"))
#   
#   results1 <- getBM(attributes = c("refsnp_id",'chr_name', "chrom_start"),
#                    filters = c("chr_name"),values = list("1"), mart = ensembl_snp)
#   # biomart query reverse strand
#    filterlist <- list(paste(credset[n,"Chr"],credset[n,"Pos"],credset[n,"Pos"],-1,sep = ":"))
#   
#   results2 <- getBM(attributes = c("refsnp_id",'chr_name', "chrom_start", "chrom_end"),
#                     filters = c("chromosomal_region"),values = filterlist, mart = ensembl_snp)
#   res=rbind(results1,results2)
#   
#   return(res)
# }
# 
# counts$Chr=gsub("\\chr*","",counts$Chr) #remove "chr" in every element in location
# 
# # load ensembl
# # ensembl = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_snp")
# # # get regions to query in correct format
# chr.region=paste(counts$Chr,counts$Start,counts$End,sep = ":")
# filterlist <- list(as.character(chr.region))
# # # biomart query
# # results <- getBM(attributes = c("refsnp_id",'chr_name', "chrom_start", "chrom_end"),
# #                  filters = c("chromosomal_region"),values = filterlist, mart = ensembl)
# 
# # biomart queries this big take too long

colnames(credset)=c("chr","start")
credset$end=credset$start

credset.gr = credset %>% 
  mutate(chr=paste0('chr', chr)) %>%   
  makeGRangesFromDataFrame

peaks.gr = makeGRangesFromDataFrame(peaks) 

# Get all credible set SNPs overlapping ATAC peaks:
overlapSNP = subsetByOverlaps(credset.gr,peaks.gr)
#overlapSNP=as.data.frame(overlapSNP)

# Get all peaks overlapping credible set SNPs :
overlapPeak = subsetByOverlaps(peaks.gr,credset.gr)
#overlapPeak=as.data.frame(overlapPeak)

# more SNPs overlaping than peaks, so there must be several SNPs within each peak in some cases

# % of SNPs from credset in peaks
(length(overlapSNP)/nrow(credset))*100

# % of peaks with intersecting SNPs from credset
(length(overlapPeak)/nrow(peaks))*100

#### Annotate by promoter, coding, etc
# 
# 
# library(ChIPpeakAnno)
# data(TSS.human.GRCh37) # load TSS info from CRCh37
# 
# # range data and annotate position with respect to closest TSS
# # This way I can see if it's in a promoter, inside a gene or in an intergenic region
# 
# range_data = function(data){
#   rangedpeak = RangedData(IRanges(start=data$Start,end=data$End), # passing peak tables as ranged data
#                           names=rownames(data),space=data$Chr)
#   return(rangedpeak)
# }
# 
# rangedPeaks = range_data(peaks)
# annotatedPeak = annotatePeakInBatch(rangedpeak, AnnotationData=TSS.human.GRCh37) # annotate with hg19. Displays distance to nearest TSS for EVERY PEAK (default: output=nearestLocation)
# annotatedPeak = annotatePeakInBatch(rangedPeaks, AnnotationData=TSS.human.GRCh37,output = "shortestDistance")
# # "shortestDistance" will output nearest features to peaks
# # Displays nearest features to TSS
# 
# annotatedPeak = as.data.frame(annotatedPeak)
source('~/WTCHG/R scripts/atac-seq/annotate_GRanges_TxDb.R') # loading whoile genome annotation function
# 
# annotate_peaks_TxDb = function(peaks.gr){
#   if(!any( is(peaks.gr)=="GRanges")){
#     stop('peak object is not GRanges')
#   }
#   
#   require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#   txdb
#   
#   # we can extract transcriptome data out of this object
#   # TX <- transcripts(txdb)
#   # 
#   # TX
#   # length(TX)
#   
#   GN = genes(txdb)
#   GN
#   length(GN)
#   
#   # intergenic
#   GN <- reduce(GN, ignore.strand=T)
#   intergenic <- gaps(GN)  # intergenic is just the complement of the genic, 
#   intergenic <- intergenic[strand(intergenic) == "*"] #This is important!!! otherwise you'll get an additional 2 entries per chromosome (one for each of + and -)
#   
#   # get promoters (certain distance from TSS)
#   PR <- promoters(txdb, upstream=2000, downstream=500)
#   PR
#   length(PR)
#   # exons
#   EX <- exons(txdb)
#   EX
#   length(EX)
#   # substract promoter parts, because they are overlapping?
#   
#   # coding - subset of exon
#   CD = cds(txdb)
#   CD
#   length(CD)
#   # substract promoter parts, because they are overlapping
#   
#   
#   # intron
#   IN = intronicParts(txdb)
#   IN
#   length(IN)
#   
#   # TX = reduce(TX, ignore.strand=T) # ignoring strand info
#   IN = reduce(IN, ignore.strand=T)
#   CD = reduce(CD, ignore.strand=T)
#   PR = reduce(PR, ignore.strand=T)
#   EX = reduce(EX, ignore.strand=T)
#   GN = reduce(GN, ignore.strand=T)
#   
#   # removing overlapping regions
#   length(intergenic)
#   intergenic = setdiff(intergenic , subsetByOverlaps( IN,intergenic))
#   length(intergenic)
#   intergenic = setdiff(intergenic , subsetByOverlaps( EX,intergenic))
#   length(intergenic)
#   intergenic = setdiff(intergenic , subsetByOverlaps( PR,intergenic))
#   length(intergenic) # slightly more sequences than before. Length!=length of sequence
#   
#   length(IN)
#   IN = setdiff(IN , subsetByOverlaps( EX,IN))
#   length(IN)
#   IN = setdiff(IN , subsetByOverlaps( PR,IN))
#   length(IN)
#   
#   length(EX)
#   EX = setdiff(EX, subsetByOverlaps(PR,EX))
#   length(EX)
#   
#   # now that I have all the features I'm interested in, I can intersect and plot the membership of 
#   # peaks
#   # peaks overlapping GWAS
#   # GWAS
#   # overall genome - compare this to published info to see I'm doing it right
#   
#   ## all peaks data
#   
#   # this way I end up with the whole peak sequence of those that overlap the features
#   PRallPeaks =subsetByOverlaps(peaks.gr,PR) # promoter
#   EXallPeaks = subsetByOverlaps(peaks.gr,EX) # exon
#   INallPeaks = subsetByOverlaps(peaks.gr,IN) # intron
#   intergenicAllPeaks = subsetByOverlaps(peaks.gr,intergenic) # intergenic
#   
#   sum(length(PRallPeaks),length(EXallPeaks),length(INallPeaks),length(intergenicAllPeaks))-
#     
#     length(subsetByOverlaps(EXallPeaks,PRallPeaks))-
#     length(subsetByOverlaps(INallPeaks,PRallPeaks))-
#     length(subsetByOverlaps(intergenicAllPeaks,PRallPeaks)) -
#     length(subsetByOverlaps(intergenicAllPeaks,EXallPeaks)) -
#     length(subsetByOverlaps(intergenicAllPeaks,INallPeaks))-
#     length(subsetByOverlaps(EXallPeaks,INallPeaks))
#   
#   # removing duplicates: prioritize assignment to PR 1st, EX 2nd, IN 3rd
#   EXallPeaks = setdiff(EXallPeaks,subsetByOverlaps(EXallPeaks,PRallPeaks)) # substracting PR peaks from the EX peaks
#   INallPeaks = setdiff(INallPeaks,subsetByOverlaps(INallPeaks,PRallPeaks)) # substracting PR peaks from the IN peaks
#   INallPeaks = setdiff(INallPeaks,subsetByOverlaps(INallPeaks,EXallPeaks)) # substracting EX peaks from the IN peaks
#   intergenicAllPeaks = setdiff(intergenicAllPeaks,subsetByOverlaps(intergenicAllPeaks,PRallPeaks)) # substracting PR peaks from the intergenic peaks
#   intergenicAllPeaks = setdiff(intergenicAllPeaks,subsetByOverlaps(intergenicAllPeaks,EXallPeaks)) # substracting EX peaks from the intergenic peaks
#   intergenicAllPeaks = setdiff(intergenicAllPeaks,subsetByOverlaps(intergenicAllPeaks,INallPeaks)) # substracting IN peaks from the intergenic peaks
#   
#   # % over 100% of peaks
#   all = sum(length(PRallPeaks),length(EXallPeaks),length(INallPeaks),length(intergenicAllPeaks))
#   
#   mylist <- list() 
#   
#   mylist[["PR"]] = (length(PRallPeaks)/all)*100
#   mylist[["EX"]] = (length(EXallPeaks)/all)*100
#   mylist[["IN"]] = (length(INallPeaks)/all)*100
#   mylist[["intergenic"]] = (length(intergenicAllPeaks)/all)*100 
# 
#   df <- do.call("rbind",mylist) #combine all vectors into a matrix
#   df = as.data.frame(df)
#   colnames(df) = "value"
#   return(df)
#   
#   
#   combined = c(PRallPeaks,EXallPeaks,INallPeaks,intergenicAllPeaks)
#   setdiff(peaks.gr,combined) # it seems like unassigned peaks are intergenic
#   
#   # safer to remove them
#   peaks.gr = setdiff(peaks.gr,setdiff(peaks.gr,combined))
#   length(peaks.gr) == length(combined)
#   
# }

# ATAC peaks and GWAS credible set SNPs
allPeaksAnnotated = annotate_GRanges_TxDb(peaks.gr) 
allPeaksAnnotated$anno = rownames(allPeaksAnnotated) # shaping for ggplot
allPeaksAnnotated$type= "All ATAC-seq peaks"
peaksOverGWASAnnotated = annotate_GRanges_TxDb(overlapPeak) 
peaksOverGWASAnnotated$anno = rownames(peaksOverGWASAnnotated)
peaksOverGWASAnnotated$type= "ATAC-seq peaks overlapping credible sets"
GWASAnnotated = annotate_GRanges_TxDb(credset.gr)
GWASAnnotated$anno = rownames(GWASAnnotated)
GWASAnnotated$type= "All SNPs in credible sets"

source('~/WTCHG/R scripts/annotate_whole_genome_TxDb.R') # loading whole genome annotation function

whole = annotate_whole_genome_TxDb() # whole genome annotation
whole$anno = rownames(whole)
whole$type = "Whole genome"
  
combined = rbind(allPeaksAnnotated,peaksOverGWASAnnotated,GWASAnnotated, whole)

combined$anno = factor(combined$anno,levels = unique(combined$anno),ordered = T,
                       labels = c("Promoter","Exon","Intron","Intergenic")) # ordering annotations
combined$type = factor(combined$type,levels = unique(combined$type),ordered = T) # ordering labels of data type

# colour scale for plot

pal=c("#96165a","#bf3b5c","#3D8DC1","#95C8EB")   # red and blue scale

source('~/WTCHG/R scripts/plotting/ggplot_stacked_barplot_percent.R') # loading barplot function

p = ggplot_stacked_barplot_percent(combined,pal,0.75 )

png("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/TxDb_annotation_peaks_credible_set_T2D_DIAGRAM.png", type="cairo",
    width=14,height=8,units="in",res=500,pointsize = 12)
p 
dev.off()


# To-do
## save peaks that overlap GWAS, and get nearest gene
## save SNPs that overlap peaks, and get nearest gene