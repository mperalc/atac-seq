# plotting longitudinal variation of peak counts


# Load necessary libraries
library(readr)  # to read big tables fast
library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(ggpubr)





# Working directory:


origin= c("Ad2.1","Ad3.1","Neo1.1")  # original samples, here called origINS

stages= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN", "BLC") # 8 differentiation stages

#stage= old_stages=c("iPSC","DE","GT","PF","PE","EN")
# origin= old_origin=c("Sbad2.1","Sbad3.4")

load("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq.xz")

cpm=as.data.frame(cpm(dge$counts))

list=c("peak_chr1_1563960","peak_chr1_975747")

############## variable: gene name and peak assigned to that gene
# list=c("INS-IGF2")
# 
# peaks=list()
# for(s in stages){
#   peaks[[s]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/DEA",
#                         "2017-12-06","_peaks_at_promoters_",s,"_across-stages.txt",sep = ""),
#                         header = T,stringsAsFactors = T)
# }
# peak_df=do.call("rbind",peaks)
# 
# peak_df_found=peak_df[match(list,peak_df$external_gene_name),]  

plot_long=cpm[match(list,rownames(cpm)),]
plot_long$Name=rownames(plot_long)
plot_long=na.omit(plot_long)              # remove NAs, in case there was not a match for a gene

diff=setdiff(list,peak_df$external_gene_name)     # which ones in the list are not in the table (probably have a different name)?


gene_number=nrow(plot_long)   # how many genes to plot

# melt data for ggplot2

long = melt(plot_long, measure.vars = c(1:(ncol(plot_long)-1)))
head(long)

# rename stages and samples
stage_2= c("iPSC", "DE", "GT", "PF", "PE", "EP","EN", "BLC")  #shortening EN names EN7= BLC (beta-like cells)
# stage_2=old_stages
samples <- c(rep(origin,8))
# samples <- c(rep(c("Sample 1","Sample 2","Sample 3"),8)) # for paper

# samples <- c(rep(origin,6))
long$variable=rep(stage_2,each=3*gene_number)                    # sample size times number of genes
# long$variable=rep(stage_2,each=2*gene_number)                    

colnames(long)[which(names(long) == "variable")] <- "stage"
long$Sample=rep(samples,each=gene_number)
long$stage <- factor(long$stage, levels=stage_2)
long$Sample=as.factor(long$Sample)
long$Name=factor(long$Name,levels=list)
head(long)


###################### plot ########################
#   
for(n in unique(long$Name)){
  
long2=long[which(long$Name==n),]
diaPalette <- c("#C15858", "#6DA567","#7883BA")  # Diabetologia palette
p <- ggplot(data = long2, aes(x = stage, y = value,group=Sample )) +
  ggtitle(unique(long2$Name)) +
  xlab ("Differentiation stages") +
  ylab ("Expression [TPM]") +
  expand_limits(y = 0) +
  geom_hline(yintercept=0,linetype="dashed",size=1,col="#DCDCDC") +
  geom_line(aes(linetype = Sample,col=Sample), size = 1) +
  scale_colour_manual(values=diaPalette) +  # diabetologia pallete
  geom_point(size=3,aes(shape=Sample,col=Sample)) +
  #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border=element_rect(size=2),axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold.italic"),
        legend.text = element_text(size=11,face="bold"),legend.title = element_text(size=13,face="bold"))

png(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/",n,".png",sep=""), type="cairo",
    width=8,height=5,units="in",res=300,pointsize = 12)
print(p)
dev.off()

}
