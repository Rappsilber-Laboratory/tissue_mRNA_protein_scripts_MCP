## Script to calculate differences in non-exponentially decaying proteins in gene pairs buffered on protein level 
## and with sustained protein co-expression
## Author: Piotr Grabowski, 2018, Rappsilber lab
library(data.table)
library(ggplot2)
library(psych)
library(reshape2)

## First get 3 lists of gene pairs:
## 1. "coreg_mRNA" = all gene pairs with co-regulated mRNAs (PCC > 0.5, BH P <0.01)
## 2. "buffered" = co-regulation buffered on protein level (mRNA PCC > 0.5, protein PCC < 0.5, mRNA BH P <0.01)
## 3. "sustained" = co-regulation sustained on protein level (mRNA PCC > 0.5, protein PCC > 0.5, BH P <0.01)
classes = fread( file.path(".","data","Gene_pair_groups_buffering.csv"))
classes[,Gene_1_sorted := ifelse(Gene_1 > Gene_2, Gene_1, Gene_2)]
classes[,Gene_2_sorted := ifelse(Gene_1 < Gene_2, Gene_1, Gene_2)] # sort row-wise
setkey(classes,Gene_1_sorted,Gene_2_sorted)

## Load NED annotation from McShane, Cell 2016
NEDs = fread( file.path(".","data","McShane_ST1_NEDs_mouse.csv"))
NEDs = NEDs[,c("Protein IDs (Uniprot)","Degradation profile")]
mapping = fread( file.path(".","data","ensembl_to_uniprot_classes_mapping.tab"))
NEDs[, ensembl := mapping$To[match(NEDs$`Protein IDs (Uniprot)`,mapping$From)]]
NEDs = NEDs[!is.na(ensembl),]
setkey(NEDs,ensembl)

## Limit our analysis to what can be annotated
classes = classes[Gene_1 %in% NEDs$ensembl & Gene_2 %in% NEDs$ensembl,]

## Add NED data to classes
classes[, Gene_1_profile := NEDs$`Degradation profile`[match(classes$Gene_1_sorted,NEDs$ensembl)]]
classes[, Gene_2_profile := NEDs$`Degradation profile`[match(classes$Gene_2_sorted,NEDs$ensembl)]]

## Keep only classes for which there are full profiles
classes_2 = classes[!is.na(Gene_1_profile) & !is.na(Gene_2_profile),]
classes_2[,min_one_NED := ifelse(Gene_1_profile == "NED" | Gene_2_profile == "NED",TRUE,FALSE)]

## Calculate NED enrichment stats
mRNA_coreg_classes = nrow(classes_2)
buffered_classes = nrow(classes_2[class == "buffered",])
sustained_classes = nrow(classes_2[class == "sustained",])

min_one_NED_counts = classes_2[,sum(min_one_NED),by=class]
sustained_counts = min_one_NED_counts[class=="sustained",V1]
buffered_counts = min_one_NED_counts[class=="buffered",V1]

## Calculate Pearson's Chi Squared test
cont_matrix = table(classes_2$class,classes_2$min_one_NED)
cont_matrix_pval = chisq.test(cont_matrix)$p.value
cont_matrix_pval

## Plot data
mRNA_coreg_NED_ratio = (sustained_counts+buffered_counts)/mRNA_coreg_classes
buffered_NED_ratio = buffered_counts/buffered_classes
sustained_NED_ratio = sustained_counts/sustained_classes

plot_data = melt(data.table(mRNA_coreg_NED_ratio,buffered_NED_ratio,sustained_NED_ratio))
plot_data$variable = factor(plot_data$variable,levels = c("mRNA_coreg_NED_ratio","buffered_NED_ratio","sustained_NED_ratio"))
plot_data = plot_data[order(plot_data$variable),]

ned_plot = ggplot(plot_data) + geom_bar(aes(x=1,y=value*100,fill=variable), stat="identity",  position = position_dodge(width=1)) + theme_bw() +
  ylab("% of pairs having at least one NED") + xlab("") +
  scale_fill_manual(values=c("steelblue","orange","darkviolet"))
ned_plot
