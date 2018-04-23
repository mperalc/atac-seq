# correlation between RNA-seq and ATAC-seq

library(readr)  # to read big tables fast
library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(ggpubr)
library(dplyr) 
library(ggpmisc)


A<-data.frame(A1=c(1,2,3,4,5),B1=c(6,7,8,9,10),C1=c(11,12,13,14,15 ))

B<-data.frame(A2=c(6,7,7,10,11),B2=c(2,1,3,8,11),C2=c(1,5,16,7,8))

cor(A,B)



origin= c("Ad2.1","Ad3.1","Neo1.1")  # original samples, here called origINS

stages= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN", "BLC") # 8 differentiation stages

# atac-seq peak counts
load("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq.xz")

cpm=as.data.frame(cpm(dge$counts)) # peak counts

#RNA counts
tpm=as.data.frame(read_tsv("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/31.01.2017.Differentiation_v2.gene.tpm.tsv"))   # read tpm file for plotting longitudinal tpm data
stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","ENstage6", "ENstage7") # 8 differentiation stages

tpm$GeneID=gsub("(ENSG[0-9]+).*", "\\1",tpm$GeneID) # cut everything after the dot


# peaks in promoters lists
peaks=list()
for(s in stages){
  peaks[[s]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/DEA",
                                   "2017-12-06","_peaks_at_promoters_",s,"_across-stages.txt",sep = ""),
                        header = T,stringsAsFactors = T)
}
peak_df=do.call("rbind",peaks)



match_DEAgenes_from_peaks_to_rna=function(peak_df,cpm,v=NA){
  # matching function: gene ID to that of RNA-seq
  # v  must be vector of gene id (ensembl)
  # otherwise takes list from atac-peak df
  if(is.na(v)){
  peak_df_found=peak_df
  }
  else{
    peak_df_found=peak_df[match(list,peak_df$ensembl_gene_id),]
    
  }
  # getting atac-seq counts
  plot_long=cpm[match(peak_df_found$Peak_names,rownames(cpm)),]
  plot_long$Peak_names=rownames(plot_long)
  plot_long=na.omit(plot_long)              # remove NAs, in case there was not a match for a gene
  plot_long=merge(peak_df_found,plot_long,by="Peak_names")
  plot_long=plot_long[c(1,12,13,18:ncol(plot_long))]
  
  ID=as.character(peak_df_found$ensembl_gene_id)
  
  pairwise_list=list() # declaring list for later
  for(p in ID){
  ## RNA matching
  
  RNA=tpm[match(p,tpm$GeneID),]  # extracts from tpm data frame the rows that contain the genes of interest
  RNA=na.omit(RNA)              # remove NAs, in case there was not a match for a gene
  if(nrow(RNA)==0){
    next
  }
  diff=setdiff(p,tpm$GeneID)           # which ones in the list are not in the table (probably have a different name)?
  
  #order the columns by stage
  
  nc=RNA[ , grepl( "iPSC" , names( RNA ) ) ]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in stage[-1])  {         
    
    i= RNA[ , grepl( s , names( RNA ) ) ]
    nc=cbind(nc,i)
  }
  
  RNA=cbind(RNA[c(1:2)],nc)
  rm(i,nc)
  
  colnames(RNA)[c(3:ncol(RNA))]=colnames(plot_long)[c(4:ncol(plot_long))]
  
  # check if there are duplicated gene ids
  if (length(which(plot_long$ensembl_gene_id==p)) > 1){
    # print(paste("There are duplicated genes from peak assignment:   ",
    #             plot_long$Peak_names, "assigned to", plot_long$external_gene_name));
    next
  }
  if (length(which(RNA$GeneID==p)) > 1){
    # print(paste("There are duplicated genes from RNA matching by Ensembl ID:    ",
    #             RNA$GeneName, "matched to", plot_long$Peak_names));
    next
  }
  unt=rbind(plot_long[which(plot_long$ensembl_gene_id==p),c(4:ncol(plot_long))],RNA[c(3:ncol(RNA))])
  rownames(unt)=c("ATAC","RNA")
  pairwise_list[[p]]=data.frame(t(unt),RNA$GeneName)
  }
  return(pairwise_list)
}

matches=match_DEAgenes_from_peaks_to_rna(peak_df,cpm)

print(paste("Found",(length(matches)/nrow(peak_df))*100,
      "% genes from the original assignation, after excluding genes that had assigned more than one peak"))

atac=as.data.frame(sapply(matches, `[[`, 1))
rna=as.data.frame(sapply(matches, `[[`, 2))
genes=sapply(matches, `[[`, 3)

correlation=diag(cor(atac,rna)) # There are missing values due to RNA-seq values that are 0 across all samples
hist(correlation)
truehist(correlation)
# p <- ggplot(pairwise_list[[ID]], aes(ATAC,RNA)) + 
#   stat_smooth(method = "lm", formula = y ~ x-1, size = 1, level=0.95) + 
#   geom_point(alpha = 0.3) +
#   stat_poly_eq(aes(label = paste(..rr.label..)),
#                label.x.npc = "right", label.y.npc = 0.15, formula = formula, 
#                parse=TRUE, size = 3) + 
#   stat_fit_glance(method = 'lm', method.args = list(formula = formula),
#                   geom='text', aes(label=ifelse(..p.value..< 0.001, "p<0.001**", 
#                                                 ifelse(..p.value..>=0.001 & ..p.value..<0.05, "p<0.05*", "p>0.05"))),
#                   label.x.npc = 'right',
#                   label.y.npc = 0.4, size = 3)
# plot(p)

png("/Users/Marta/Documents/WTCHG/DPhil/Plots/INS-IGF2_atac_rna_cor.png", type="cairo", antialias="default",
    width=5,height=5,units="in",res=600,pointsize = 13)

corrplot=ggscatter(pairwise_list[[ID]], x = "ATAC", y = "RNA", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "ATAC-seq counts (cpm) for closest peak in promoter ",
                   ylab = "Gene expression (TPM)",
                   title = as.character(pairwise_list[[ID]]$list[[1]]))
plot(corrplot)
dev.off()
#### other possible gene
name="INS"

RNA=tpm[match(name,tpm$GeneName),]  # extracts from tpm data frame the rows that contain the genes of interest
RNA=na.omit(RNA)              # remove NAs, in case there was not a match for a gene

#order the columns by stage

nc=RNA[ , grepl( "iPSC" , names( RNA ) ) ]  #takes columns whose name matches x

#Do the same for the other stages

for (s in stage[-1])  {         
  
  i= RNA[ , grepl( s , names( RNA ) ) ]
  nc=cbind(nc,i)
}

RNA=cbind(RNA[c(1:2)],nc)
rm(i,nc)


colnames(RNA)[c(3:ncol(RNA))]=colnames(plot_long)[c(4:ncol(plot_long))]
unt=rbind(plot_long[c(4:ncol(plot_long))],RNA[c(3:ncol(RNA))])
rownames(unt)=c("ATAC","RNA")
pairwise_list=list()
pairwise_list[[name]]=data.frame(t(unt),name)

png("/Users/Marta/Documents/WTCHG/DPhil/Plots/INS_atac_rna_cor.png", type="cairo", antialias="default",
    width=5,height=5,units="in",res=600,pointsize = 13)

corrplot=ggscatter(pairwise_list[[name]], x = "ATAC", y = "RNA", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "ATAC-seq counts (cpm) for closest peak in promoter ",
                   ylab = "Gene expression (TPM)",
                   title = as.character(pairwise_list[[name]]$name[[1]]))


plot(corrplot)
dev.off()
