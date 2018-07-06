## Script to calculate differences in transcription and translation rates from Selbach 2011 in gene pairs buffered on protein level 
## and with sustained protein co-expression
## Author: Piotr Grabowski, 2018, Rappsilber lab
library(data.table)
library(ggplot2)
library(gridExtra)
library(psych)
library(reshape2)
library(biomaRt)
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl) ## sometimes one has to retry multiple times to load it...

## First get 3 lists of gene pairs:
## 1. "coreg_mRNA" = all gene pairs with co-regulated mRNAs (PCC > 0.5, BH P <0.01)
## 2. "buffered" = co-regulation buffered on protein level (mRNA PCC > 0.5, protein PCC < 0.5, mRNA BH P <0.01)
## 3. "sustained" = co-regulation sustained on protein level (mRNA PCC > 0.5, protein PCC > 0.5, BH P <0.01)
classes = fread( file.path(".","data","Gene_pair_groups_buffering.csv"))
classes[,Gene_1_sorted := ifelse(Gene_1 > Gene_2, Gene_1, Gene_2)]
classes[,Gene_2_sorted := ifelse(Gene_1 < Gene_2, Gene_1, Gene_2)] # sort row-wise
setkey(classes,Gene_1_sorted,Gene_2_sorted)

## Load Schwannhaeusser 2011 Nature data
rates = fread( file.path(".","data","schwannhaeusser_selbach_2011_rates.csv"))
rates[,first_ensembl_tx := tstrsplit(`ENSEMBL ID`,";",keep=1)]
mapping = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id"),filters = "ensembl_transcript_id", values = rates$first_ensembl_tx, mart=ensembl)
rates[,ensembl_gene_id := mapping$ensembl_gene_id[match(rates$first_ensembl_tx,mapping$ensembl_transcript_id)]]
rates = rates[,.(ensembl_gene_id,`transcription rate (vsr) average [molecules/(cell*h)]`,`translation rate constant (ksp) average [molecules/(mRNA*h)]` )]
rates = rates[!is.na(ensembl_gene_id),]

## Add the translation rates
classes[,Gene_1_tl := rates$`translation rate constant (ksp) average [molecules/(mRNA*h)]`[match(Gene_1_sorted,rates$ensembl_gene_id)]]
classes[,Gene_2_tl := rates$`translation rate constant (ksp) average [molecules/(mRNA*h)]`[match(Gene_2_sorted,rates$ensembl_gene_id)]]


############## translation check ############## 
classes_tl = classes[!is.na(Gene_1_tl) & !is.na(Gene_2_tl),]
classes_tl[,tl_ratio := abs(log2(Gene_1_tl/Gene_2_tl))]
## Tag a pair with having similar tl ratio if the value is < 1
classes_tl[,sim_tl_rate := ifelse(tl_ratio <= 1, TRUE, FALSE)]
classes_tl_cont = table(classes_tl$class,classes_tl$sim_tl_rate)
tl_pval = chisq.test(classes_tl_cont)$p.value
tl_plot_ratios = classes_tl_cont[,2]/(classes_tl_cont[,1]+classes_tl_cont[,2])
tl_plot_ratios = data.frame(tl_plot_ratios)
tl_plot_ratios$class = row.names(tl_plot_ratios)
tl_plot_ratios["mRNA_coreg",] = c(sum(classes_tl_cont[,2]) / (sum(classes_tl_cont[,2]) + sum(classes_tl_cont[,1])),"mRNA_coreg")
tl_plot_ratios$numbers = c(nrow(classes_tl[class == "buffered",]),
                           nrow(classes_tl[class == "sustained",]),nrow(classes_tl))
tl_plot_ratios$tl_plot_ratios = as.numeric(tl_plot_ratios$tl_plot_ratios)

## order
tl_plot_ratios$class = factor(tl_plot_ratios$class,levels=c("mRNA_coreg","buffered","sustained"))
tl_plot_ratios = tl_plot_ratios[order(tl_plot_ratios$class),]
tl_plot_ratios

## Plot data
sim_tl_rates_plot = ggplot(tl_plot_ratios) + geom_bar(aes(x=1,y=tl_plot_ratios*100,fill=class), stat="identity",  position = position_dodge(width=1)) + theme_bw() +
  ylab("% of pairs having similar translation rates") + xlab("") +
  scale_fill_manual(values=c("steelblue","orange","darkviolet"))
sim_tl_rates_plot
