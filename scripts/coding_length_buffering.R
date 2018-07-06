## Script to calculate differences in CDS lenghts in gene pairs buffered on protein level 
## and with sustained protein co-expression
## Author: Piotr Grabowski, 2018, Rappsilber lab
library(data.table)
library(ggplot2)
library(seqinr)
library(reshape2)
library(plyr)

## First get 3 lists of gene pairs:
## 1. "coreg_mRNA" = all gene pairs with co-regulated mRNAs (PCC > 0.5, BH P <0.01)
## 2. "buffered" = co-regulation buffered on protein level (mRNA PCC > 0.5, protein PCC < 0.5, mRNA BH P <0.01)
## 3. "sustained" = co-regulation sustained on protein level (mRNA PCC > 0.5, protein PCC > 0.5, BH P <0.01)
classes = fread( file.path(".","data","Gene_pair_groups_buffering.csv"))
classes[,Gene_1_sorted := ifelse(Gene_1 > Gene_2, Gene_1, Gene_2)]
classes[,Gene_2_sorted := ifelse(Gene_1 < Gene_2, Gene_1, Gene_2)] # sort row-wise
setkey(classes,Gene_1_sorted,Gene_2_sorted)

## Load CDS sequences and calculate lengths
cds = read.fasta( file.path(".","data","mouse_coding_seqs_biomart_25052018.txt"))
cds = ldply(lapply(cds,length))

## Limit classes data to whatever we have annotation for
classes = classes[Gene_1 %in% cds[,1] & Gene_2 %in% cds[,1],]

## Coding length is considered similar if the longer protein was < 1.5‐fold longer than the shorter protein. 
## Add cds lengths to classes
classes$Gene_1_sorted_len = cds$V1[match(classes$Gene_1_sorted,cds$.id)]
classes$Gene_2_sorted_len = cds$V1[match(classes$Gene_2_sorted,cds$.id)]
classes[, cds_len_ratio := ifelse(Gene_1_sorted_len > Gene_2_sorted_len, 
                                  Gene_1_sorted_len/Gene_2_sorted_len,
                                  Gene_2_sorted_len/Gene_1_sorted_len)]
classes[, similar_len := ifelse(cds_len_ratio < 1.5, TRUE, FALSE)]

## Calculate similar lenghts enrichment stats
mRNA_coreg_classes = nrow(classes)
buffered_classes = nrow(classes[class == "buffered",])
sustained_classes = nrow(classes[class == "sustained",])
similar_len_counts = classes[,sum(similar_len,na.rm=T),by=class]
sustained_counts = similar_len_counts[class=="sustained",V1]
buffered_counts = similar_len_counts[class=="buffered",V1]

## Calculate Pearson's Chi Squared
cont_matrix = table(classes$similar_len,classes$class)
cont_matrix_pval = chisq.test(cont_matrix)$p.value

## Plot data
mRNA_coreg_similar_len_ratio = (sustained_counts+buffered_counts)/mRNA_coreg_classes
buffered_similar_len_ratio = buffered_counts/buffered_classes
sustained_similar_len_ratio = sustained_counts/sustained_classes

plot_data = melt(data.table(mRNA_coreg_similar_len_ratio,buffered_similar_len_ratio,sustained_similar_len_ratio))
cod_len_plot = ggplot(plot_data) + geom_bar(aes(x=1,y=value*100,fill=variable), stat="identity",  position = position_dodge(width=1)) + theme_bw() +
  ylab("% of pairs having similar coding lengths") + xlab("") +
  scale_fill_manual(values=c("steelblue","orange","darkviolet"))
cod_len_plot

