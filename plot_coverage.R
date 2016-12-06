# plot Coverage from atac-seq bam files

#  that is, average overlapping reads per bp


setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/coverage")


currentDate <- Sys.Date() # to save date in name of output files


table <- read.table(file="coverage.tab") # tab separated file with coverage from 1 million regions sampled randomly

colnames(table)[1:3]=c("chr","start","end")


samples=c("A","B","C","D")
samples_2=rep(samples,each=6)
stage=as.character(c(1:6))
stage_2=rep(stage,4)
colnames(table)[4:ncol(table)]=paste(samples_2,stage_2,sep="")   # samples created

coverage=data.frame(matrix(nrow=24,ncol=2))
coverage$X1=paste(samples_2,stage_2,sep="")
coverage$X2=as.numeric(colMeans(table[4:ncol(table)]))   # table with means and samples created

colnames(coverage)=c("samples","mean_cov")
coverage$mean_cov=round(coverage$mean_cov,2)   # rounding

# plot the histogram of number of reads aligned to each bp

## do

p <- ggplot(data = long, aes(x = stage, y = value,group=sample )) +
  ggtitle("Counts per gene, stage & sample") +
  xlab ("Differentiation stages") +
  ylab ("Transcripts per million") +
  geom_line(aes(col = sample), size = 1) + facet_wrap(~GeneName,scales="free")
print(p)

