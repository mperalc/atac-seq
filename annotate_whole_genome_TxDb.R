# Annotate whole genome with TxDb
# get percentage of genome that belongs to each annotation
# promoter, intron, exon, intergenic
annotate_whole_genome_TxDb = function(){
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  # GN = genes(txdb)
  # GN
  # length(GN)
  # 
  # # intergenic
  # #GN <- reduce(GN, ignore.strand=T)
  # intergenic <- gaps(GN)  # intergenic is just the complement of the genic, 
  # intergenic <- intergenic[strand(intergenic) == "+" | strand(intergenic) == "-" ] #
  # 
  # # get promoters (certain distance from TSS)
  # PR <- promoters(txdb, upstream=2000, downstream=500)
  # PR
  # length(PR)
  # # exons
  # EX <- exons(txdb)
  # EX
  # length(EX)
  # # substract promoter parts, because they are overlapping?
  # 
  # # coding - subset of exon
  # CD = cds(txdb)
  # CD
  # length(CD)
  # # substract promoter parts, because they are overlapping
  # 
  # # intron
  # IN = intronicParts(txdb)
  # IN
  # length(IN)
  # 
  # length(intergenic)
  # intergenic = setdiff(intergenic , subsetByOverlaps( IN,intergenic))
  # length(intergenic)
  # intergenic = setdiff(intergenic , subsetByOverlaps( EX,intergenic))
  # length(intergenic)
  # intergenic = setdiff(intergenic , subsetByOverlaps( PR,intergenic))
  # length(intergenic) # slightly more sequences than before. Length!=length of sequence
  # 
  # length(IN)
  # IN = setdiff(IN , subsetByOverlaps( EX,IN))
  # length(IN)
  # IN = setdiff(IN , subsetByOverlaps( PR,IN))
  # length(IN)
  # 
  # length(EX)
  # EX = setdiff(EX, subsetByOverlaps(PR,EX))
  # length(EX)
  # 
  # ######################### whole genome or just autosomes?
  # EX_sum = sum(as.numeric(width(EX)))
  # PR_sum = sum(as.numeric(width(PR)))
  # IN_sum = sum(as.numeric(width(IN)))
  # intergenic_sum = sum(as.numeric(width(intergenic)))
  # total_sum = EX_sum + PR_sum + IN_sum + intergenic_sum  # this includes *both* strands, so half of that would be the length of genome
  # (total_sum / 2)/1000000
  # # my estimate is a bit under the latest reports (3,234.83 Mb)
  # 
  # # In this case the % of each annotation has been calculated by total bp and not by feature count (ATAC, SNP... binary measure not taking into account their size)
  # # because obviously the intergenic regions for example are fewer than exons, but span more basepairs and therefore can potentially overlap more features
  # # I believe that coding the ATAC by overlap w/o taking into account their length and assigning only to one feature 
  # # (potentially the one where it would have more effect) facilitates the interpretation of the data and reflects the biology better
  # # eg peaks overlapping promoters and introns are more likely to act on promoters, as open access for TFs
  # 
  # mylist <- list() 
  # 
  # mylist[["PR"]] = (PR_sum/total_sum)*100
  # mylist[["EX"]] = (EX_sum/total_sum)*100
  # mylist[["IN"]] = (IN_sum/total_sum)*100
  # mylist[["intergenic"]] = (intergenic_sum/total_sum)*100
  # 
  # df <- do.call("rbind",mylist) #combine all vectors into a matrix
  # df = as.data.frame(df)
  # colnames(df) = "value"
  
  # results here give very low intron and very high intergenic
  # unstranded calculations more similar to published data
  # could the unstranded method be overestimating the amount of introns by counting complementary strands where the technical definition should be intergenic?
  # Maybe the unstranded definition of annotations is more biologically relevant?
  
  # return(df)
  
  ############3############ unstranded  #################################
  GN = genes(txdb)
  GN
  length(GN)
  
  # intergenic
  GN <- reduce(GN, ignore.strand=T)
  intergenic <- gaps(GN)  # intergenic is just the complement of the genic, 
  intergenic <- intergenic[strand(intergenic) == "*"] #This is important!!! otherwise you'll get an additional 2 entries per chromosome (one for each of + and -)
  
  # get promoters (certain distance from TSS)
  PR <- promoters(txdb, upstream=2000, downstream=500)
  PR
  length(PR)
  # exons
  EX <- exons(txdb)
  EX
  length(EX)
  # substract promoter parts, because they are overlapping?
  
  # coding - subset of exon
  CD = cds(txdb)
  CD
  length(CD)
  # substract promoter parts, because they are overlapping
  
  
  # intron
  IN = intronicParts(txdb)
  IN
  length(IN)
  
  # TX = reduce(TX, ignore.strand=T) # ignoring strand info
  IN = reduce(IN, ignore.strand=T)
  CD = reduce(CD, ignore.strand=T)
  PR = reduce(PR, ignore.strand=T)
  EX = reduce(EX, ignore.strand=T)
  GN = reduce(GN, ignore.strand=T)
  
  # removing overlapping regions
  length(intergenic)
  intergenic = setdiff(intergenic , subsetByOverlaps( IN,intergenic))
  length(intergenic)
  intergenic = setdiff(intergenic , subsetByOverlaps( EX,intergenic))
  length(intergenic)
  intergenic = setdiff(intergenic , subsetByOverlaps( PR,intergenic))
  length(intergenic) # slightly more sequences than before. Length!=length of sequence
  
  length(IN)
  IN = setdiff(IN , subsetByOverlaps( EX,IN))
  length(IN)
  IN = setdiff(IN , subsetByOverlaps( PR,IN))
  length(IN)
  
  length(EX)
  EX = setdiff(EX, subsetByOverlaps(PR,EX))
  length(EX)
  
  EX_sum = sum(as.numeric(width(EX)))
  PR_sum = sum(as.numeric(width(PR)))
  IN_sum = sum(as.numeric(width(IN)))
  intergenic_sum = sum(as.numeric(width(intergenic)))
  total_sum = EX_sum + PR_sum + IN_sum + intergenic_sum 
  # my estimate is a bit under the latest reports (3,234.83 Mb)
  
  # In this case the % of each annotation has been calculated by total bp and not by feature count (ATAC, SNP... binary measure not taking into account their size)
  # because obviously the intergenic regions for example are fewer than exons, but span more basepairs and therefore can potentially overlap more features
  # I believe that coding the ATAC by overlap w/o taking into account their length and assigning only to one feature 
  # (potentially the one where it would have more effect) facilitates the interpretation of the data and reflects the biology better
  # eg peaks overlapping promoters and introns are more likely to act on promoters, as open access for TFs
  
  # for this reason promoters could be overestimated here (vs estimate of 0.5% by Alexander et al., Nat Rev Genet, 2010)
  
  mylist <- list() 
  
  mylist[["PR"]] = (PR_sum/total_sum)*100
  mylist[["EX"]] = (EX_sum/total_sum)*100
  mylist[["IN"]] = (IN_sum/total_sum)*100
  mylist[["intergenic"]] = (intergenic_sum/total_sum)*100
  
  df <- do.call("rbind",mylist) #combine all vectors into a matrix
  df = as.data.frame(df)
  colnames(df) = "value"
  return(df)
}
