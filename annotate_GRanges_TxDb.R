# takes in object in GRanges format and returns dataframe with percentages of promoters, exons, introns and intergenic
annotate_GRanges_TxDb = function(gr){
  if(!any( is(gr)=="GRanges")){
    stop('peak object is not GRanges')
  }
  
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  txdb
  
  # we can extract transcriptome data out of this object
  # TX <- transcripts(txdb)
  # 
  # TX
  # length(TX)
  
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
  
  # now that I have all the features I'm interested in, I can intersect and plot the membership of 
  # peaks
  # peaks overlapping GWAS
  # GWAS
  # overall genome - compare this to published info to see I'm doing it right
  
  ## all peaks data
  
  # this way I end up with the whole peak sequence of those that overlap the features
  PRallPeaks =subsetByOverlaps(gr,PR) # promoter
  EXallPeaks = subsetByOverlaps(gr,EX) # exon
  INallPeaks = subsetByOverlaps(gr,IN) # intron
  intergenicAllPeaks = subsetByOverlaps(gr,intergenic) # intergenic
  
  sum(length(PRallPeaks),length(EXallPeaks),length(INallPeaks),length(intergenicAllPeaks))-
    
    length(subsetByOverlaps(EXallPeaks,PRallPeaks))-
    length(subsetByOverlaps(INallPeaks,PRallPeaks))-
    length(subsetByOverlaps(intergenicAllPeaks,PRallPeaks)) -
    length(subsetByOverlaps(intergenicAllPeaks,EXallPeaks)) -
    length(subsetByOverlaps(intergenicAllPeaks,INallPeaks))-
    length(subsetByOverlaps(EXallPeaks,INallPeaks))
  
  # removing duplicates: prioritize assignment to PR 1st, EX 2nd, IN 3rd
  EXallPeaks = setdiff(EXallPeaks,subsetByOverlaps(EXallPeaks,PRallPeaks)) # substracting PR peaks from the EX peaks
  INallPeaks = setdiff(INallPeaks,subsetByOverlaps(INallPeaks,PRallPeaks)) # substracting PR peaks from the IN peaks
  INallPeaks = setdiff(INallPeaks,subsetByOverlaps(INallPeaks,EXallPeaks)) # substracting EX peaks from the IN peaks
  intergenicAllPeaks = setdiff(intergenicAllPeaks,subsetByOverlaps(intergenicAllPeaks,PRallPeaks)) # substracting PR peaks from the intergenic peaks
  intergenicAllPeaks = setdiff(intergenicAllPeaks,subsetByOverlaps(intergenicAllPeaks,EXallPeaks)) # substracting EX peaks from the intergenic peaks
  intergenicAllPeaks = setdiff(intergenicAllPeaks,subsetByOverlaps(intergenicAllPeaks,INallPeaks)) # substracting IN peaks from the intergenic peaks
  
  # % over 100% of peaks
  all = sum(length(PRallPeaks),length(EXallPeaks),length(INallPeaks),length(intergenicAllPeaks))
  
  mylist <- list() 
  
  mylist[["PR"]] = (length(PRallPeaks)/all)*100
  mylist[["EX"]] = (length(EXallPeaks)/all)*100
  mylist[["IN"]] = (length(INallPeaks)/all)*100
  mylist[["intergenic"]] = (length(intergenicAllPeaks)/all)*100 
  
  df <- do.call("rbind",mylist) #combine all vectors into a matrix
  df = as.data.frame(df)
  colnames(df) = "value"
  return(df)
  
  
  combined = c(PRallPeaks,EXallPeaks,INallPeaks,intergenicAllPeaks)
  setdiff(gr,combined) # it seems like unassigned peaks are intergenic
  
  # safer to remove them
  gr = setdiff(gr,setdiff(gr,combined))
  length(gr) == length(combined)
  
}
