## Script to calculate differences in ribosome profiles in gene pairs buffered on protein level and with sustained protein co-expression
## Author: Piotr Grabowski, 2018, Rappsilber lab
library(data.table)
library(ggplot2)
library(reshape2)

## First get 3 lists of gene pairs:
## 1. "coreg_mRNA" = all gene pairs with co-regulated mRNAs (PCC > 0.5, BH P <0.01)
## 2. "buffered" = co-regulation buffered on protein level (mRNA PCC > 0.5, protein PCC < 0.5, mRNA BH P <0.01)
## 3. "sustained" = co-regulation sustained on protein level (mRNA PCC > 0.5, protein PCC > 0.5, BH P <0.01)
classes = fread( file.path(".","data","Gene_pair_groups_buffering.csv"))
classes[,Gene_1_sorted := ifelse(Gene_1 > Gene_2, Gene_1, Gene_2)]
classes[,Gene_2_sorted := ifelse(Gene_1 < Gene_2, Gene_1, Gene_2)] # sort row-wise
setkey(classes,Gene_1_sorted,Gene_2_sorted)

## Load a list of gene pairs available in the ribosome profiling exp
genes_Gatfield = fread( file.path(".","data","Gatfield_liver_processed_genes_available.csv"), header=F)$V1

## Load a list of gene pairs with correlated ribosome profiles
riboprofiles = fread( file.path(".","data","Gatfield_liver_processed_correlated.csv"))
riboprofiles[,Gene_1_sorted := ifelse(Gene_1 > Gene_2, Gene_1, Gene_2)]
riboprofiles[,Gene_2_sorted := ifelse(Gene_1 < Gene_2, Gene_1, Gene_2)] # sort row-wise
setkey(riboprofiles,Gene_1_sorted,Gene_2_sorted)

## Limit classes list to whatever we have annotation for
genes_with_ribo_profiles = unique( c(riboprofiles$Gene_1,riboprofiles$Gene_2) )
classes = classes[Gene_1 %in% genes_with_ribo_profiles & Gene_2 %in% genes_with_ribo_profiles,]

## Merge data
merged = riboprofiles[classes]
merged$rp_correlated = FALSE
merged$rp_correlated[!is.na(merged$PCC_ribo_profiles)] = TRUE

## Calculate counts for stat analysis and plotting
mRNA_coreg_classes = nrow(classes)
mRNA_coreg_ratio = nrow(merged[rp_correlated==TRUE & class != "prot_coreg_only"])/nrow(merged)
buffered_ratio = nrow(merged[rp_correlated==TRUE & class=="buffered",])/nrow(merged[class=="buffered",])
sust_ratio = nrow(merged[rp_correlated==TRUE & class=="sustained",])/nrow(merged[class=="sustained",])

## 3x2 Contigency for Pearson's Chi Squared test
buf_vs_sust = table(merged$class,merged$rp_correlated)
buf_vs_sust_chi = chisq.test(buf_vs_sust)$p.value

## Plot
ratios  = c(mRNA_coreg_ratio,buffered_ratio,sust_ratio)
gene_pair_numbers = c(nrow(merged),nrow(merged[class=="buffered",]),nrow(merged[class=="sustained",]))
plot_data = data.frame(class = c("mRNA_coreg","buffered","sustained"),ratios,gene_pair_numbers)
plot_data$class = factor(c("mRNA_coreg","buffered","sustained"), levels = c("mRNA_coreg","buffered","sustained"))

ribo_profile_plot = ggplot(plot_data) + geom_bar(aes(x=1,y=ratios*100,fill=class),stat = "identity", position = position_dodge(width=1)) + 
  theme_bw() + scale_fill_manual(values=c("steelblue","orange","darkviolet")) +
  ylab("% of pairs with correlated ribosome profiles") +
  xlab("") + ylim(0,12)
ribo_profile_plot
