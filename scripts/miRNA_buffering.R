## Script to calculate differences in miRNA targeting in gene pairs buffered on protein level 
## and with sustained protein co-expression
## Author: Piotr Grabowski, 2018, Rappsilber lab
library(data.table)
library(ggplot2)
library(psych)
library(reshape2)
library(biomaRt)
library(cowplot)
library(plyr)
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)

## First get 3 lists of gene pairs:
## 1. "coreg_mRNA" = all gene pairs with co-regulated mRNAs (PCC > 0.5, BH P <0.01)
## 2. "buffered" = co-regulation buffered on protein level (mRNA PCC > 0.5, protein PCC < 0.5, mRNA BH P <0.01)
## 3. "sustained" = co-regulation sustained on protein level (mRNA PCC > 0.5, protein PCC > 0.5, BH P <0.01)
classes = fread( file.path(".","data","Gene_pair_groups_buffering.csv"))
classes[,Gene_1_sorted := ifelse(Gene_1 > Gene_2, Gene_1, Gene_2)]
classes[,Gene_2_sorted := ifelse(Gene_1 < Gene_2, Gene_1, Gene_2)] # sort row-wise
setkey(classes,Gene_1_sorted,Gene_2_sorted)

## 1. Load mouse miRNA data and check if proteins with sustained protein coregulation tend to be targets of the same miRNAs
mirna = fread(file.path(".","data","mouse_brain_CLEAR_miRNAs_Darnell2015.csv"))
mapping_table = getBM(attributes = c("entrezgene","ensembl_gene_id"),filters = "entrezgene",values=unique(mirna$gene.id), mart=ensembl)
mirna[,ensembl := mapping_table$ensembl_gene_id[match(gene.id,mapping_table$entrezgene)]]
mirna = mirna[!is.na(ensembl),]
mirna_full = copy(mirna) ## make a copy for later to avoid reloading
mirna = mirna[grep("UTR|CDS",region),]
mirna = mirna[grep("CDS",region),] ## subset only to CDS
mirna = mirna[-grep("UTR",region),] ## subset only to CDS
mirna = unique(mirna[,.(miRNA,ensembl)]) 

## create a data table of genes sharing miRNAs
interactions = data.table(dplyr::inner_join(mirna,mirna,by="miRNA"))
interactions = interactions[ensembl.x != ensembl.y, ]

## Limit classes list to whatever we have annotation for
classes = classes[Gene_1 %in% mirna$ensembl & Gene_2 %in% mirna$ensembl,]

## sort ensembl ids row-wise
interactions[,Gene_1_sorted := ifelse(ensembl.x > ensembl.y, ensembl.x, ensembl.y)]
interactions[,Gene_2_sorted := ifelse(ensembl.x < ensembl.y, ensembl.x, ensembl.y)]
interactions = interactions[, .(miRNA, Gene_1_sorted, Gene_2_sorted)]
interactions = unique(interactions)
setkey(interactions,Gene_1_sorted,Gene_2_sorted)

## Calculate how many common miRNAs per gene pair
gene_pair_mirna_counts = data.table(plyr::count(interactions[,.(Gene_1_sorted,Gene_2_sorted)]))
gene_pair_mirna_counts = gene_pair_mirna_counts[order(gene_pair_mirna_counts$freq, decreasing = T),]
setkey(gene_pair_mirna_counts, Gene_1_sorted,Gene_2_sorted)

## Define groups as targeted by the same miRNA if there is at least one miRNA which targets both
shared_targets = merge(gene_pair_mirna_counts,classes,all.y=T)
shared_targets = shared_targets[!duplicated(shared_targets[,.(Gene_1,Gene_2)]),]
shared_targets$freq[is.na(shared_targets$freq)] = 0

## Look at the mean numbers of shared miRNAs per class
mean_number_shared = shared_targets[,mean(freq),by=class]
colnames(mean_number_shared)[2] = "mean_shared_miRNAs"
mean_number_shared = rbind( mean_number_shared,list("mrna_coreg",mean(shared_targets$freq)) )
mean_number_shared$class = factor(mean_number_shared$class,levels=c("mrna_coreg","buffered","sustained"))
mean_number_shared = mean_number_shared[order(class),]

## check how significant is that mean shift
mann_whitney_pval = wilcox.test(shared_targets[class=="buffered",freq],shared_targets[class=="sustained",freq])$p.value
mann_whitney_pval

cds_buff_plot = ggplot(mean_number_shared) + geom_bar(aes(x=1,y=mean_shared_miRNAs,fill=class), stat="identity",  position = position_dodge(width=1)) +
  theme_bw() +
  scale_fill_manual(values=c("steelblue","orange","darkviolet")) +
  ylab("Mean number of shared miRNAs per gene pair") + theme(legend.position = "none") + xlab("")
cds_buff_plot

