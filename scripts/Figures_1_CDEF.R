## Script for generating Figures 1 C,D,E,F in the manuscript
## Authors: Piotr Grabowski and Georg Kustatscher, 05.2018

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)

setwd("") ## set the working directory, the data will be loaded using relative paths

############# Figure 1C consecutive reactions #############
# Load data
df = read.csv( file.path("data","mouse_SILAC_TPMs_log2_final_min8_features.csv"), stringsAsFactors=FALSE)
df = df[ !is.na( df[, "Ensembl"] ) ,]
rownames(df) = df[, "Ensembl"]
TPM = df[, grep("TPM_", colnames(df))]
SILAC = df[, grep("SILAC_", colnames(df))]
tTPM = t(TPM)
TPM_tissue_medians = apply(tTPM, 1, median, na.rm=TRUE)
tTPM_mn = sweep(tTPM, 1, TPM_tissue_medians, FUN="-")
tSILAC = t(SILAC)
SILAC_tissue_medians = apply(tSILAC, 1, median, na.rm=TRUE)
tSILAC_mn = sweep(tSILAC, 1, SILAC_tissue_medians, FUN="-")

# Load test data
all_test_pairs = read.csv( file.path("data","Reactome_consec_reacts_mouse.csv"), stringsAsFactors = FALSE)
is.na(all_test_pairs$Gene_1) = !(all_test_pairs$Gene_1 %in% rownames(df)) # turn NA if gene is NOT in our data
is.na(all_test_pairs$Gene_2) = !(all_test_pairs$Gene_2 %in% rownames(df)) # turn NA if gene is NOT in our data
our_gene_pairs = all_test_pairs[complete.cases(all_test_pairs),]   # remove rows with NA

# Compile matrices with input data for Gene_1 and Gene_2, respectively
TPM_Gene_1 = tTPM_mn[, match(our_gene_pairs$Gene_1, colnames(tTPM_mn))]
TPM_Gene_2 = tTPM_mn[, match(our_gene_pairs$Gene_2, colnames(tTPM_mn))]
SILAC_Gene_1 = tSILAC_mn[, match(our_gene_pairs$Gene_1, colnames(tSILAC_mn))]
SILAC_Gene_2 = tSILAC_mn[, match(our_gene_pairs$Gene_2, colnames(tSILAC_mn))]
TPM_cor = numeric()
SILAC_cor = numeric()

for (i in 1:nrow(our_gene_pairs)) {
  temp_TPM_1 = TPM_Gene_1[,i]
  temp_TPM_2 = TPM_Gene_2[,i]
  TPM_cor[i] = cor(temp_TPM_1,temp_TPM_2, use= "pairwise.complete.obs",  method= c("pearson"))
  temp_SILAC_1 = SILAC_Gene_1[,i]
  temp_SILAC_2 = SILAC_Gene_2[,i]
  SILAC_cor[i] = cor(temp_SILAC_1,temp_SILAC_2, use= "pairwise.complete.obs",  method= c("pearson"))
}
gene_pairs_cor = data.frame(Gene_1 = colnames(TPM_Gene_1), Gene_2 = colnames(TPM_Gene_2),
                             mRNA_correlation = TPM_cor, Protein_correlation = SILAC_cor)  
gene_pairs_cor = merge(gene_pairs_cor, our_gene_pairs, sort=FALSE)

# Create control by shuffling gene pairs

# Conditional shuffling here means that the test genes are paired randomly, but under the following conditions:
# (a) genes are never randomly paired with themselves
# (b) all randomly generated gene pairs are unique
# IMPORTANT: Re-execute this section if you receive an error message. The error occurs only if no more random pairs
# can be created due to the above restrictions.
gene_pool = c(our_gene_pairs$Gene_1, our_gene_pairs$Gene_2)             # Create a pool with all test genes

random_pairs_1 = character()                                            # Create empty vectors to fill up in loop
random_pairs_2 = character() 

while ( length(gene_pool) > 0) {
  
  Random_Gene_1 = sample( gene_pool ,1)                                   # Randomly choose a new Gene_1 
  
  previous_match_A = random_pairs_2[ random_pairs_1 == Random_Gene_1 ]  # When new Gene_1 was chosen as Random_gene_1 in previous iteration, which genes were paired as Random_gene_2?
  previous_match_B = random_pairs_1[ random_pairs_2 == Random_Gene_1 ]  # When new Gene_1 was chosen as Random_gene_2 in previous iteration, which genes were paired as Random_gene_1?
  previous_matches = c(previous_match_A, previous_match_B)              # All genes that were randomly matched to the current gene 1 in previous iterations
  previous_matches_indices = which( gene_pool %in% previous_matches )   # Indices of all genes that were previously matched to current gene 1
  
  self_indices = which( gene_pool == Random_Gene_1 )                    # Indices of all occurrences of Random_Gene_1, to avoid pairing it up with itself
  
  exclude = c(previous_matches_indices, self_indices)                   # Genes that must not be drawn as random gene 2
  
  Random_Gene_2 = sample( gene_pool[-exclude] ,1)                       # Randomly choose a new Gene_2 (with exceptions)
  
  random_pairs_1 = c(random_pairs_1, Random_Gene_1)                     # Add the randomly chosen genes to the appropriate vectors
  random_pairs_2 = c(random_pairs_2, Random_Gene_2)
  
  gene_pool = gene_pool[ -match( c(Random_Gene_1, Random_Gene_2) , gene_pool) ]    # Remove randomly chosen gene pair from the pool
}

random_gene_pairs = data.frame(Random_Gene_1=random_pairs_1, Random_Gene_2=random_pairs_2, stringsAsFactors = FALSE)    # Create dataframe with random gene pairs

# Compile matrices with input data for Random_Gene_1 and Random_Gene_2, respectively
random_TPM_Gene_1 = tTPM_mn[, match(random_gene_pairs$Random_Gene_1, colnames(tTPM_mn))]
random_TPM_Gene_2 = tTPM_mn[, match(random_gene_pairs$Random_Gene_2, colnames(tTPM_mn))]
random_SILAC_Gene_1 = tSILAC_mn[, match(random_gene_pairs$Random_Gene_1, colnames(tSILAC_mn))]
random_SILAC_Gene_2 = tSILAC_mn[, match(random_gene_pairs$Random_Gene_2, colnames(tSILAC_mn))]
random_TPM_cor = numeric()
random_SILAC_cor = numeric()

for (i in 1:nrow(random_gene_pairs)) {
  temp_random_TPM_1 = random_TPM_Gene_1[,i]
  temp_random_TPM_2 = random_TPM_Gene_2[,i]
  random_TPM_cor[i] = cor(temp_random_TPM_1,temp_random_TPM_2, use= "pairwise.complete.obs",  method= c("pearson"))
  temp_random_SILAC_1 = random_SILAC_Gene_1[,i]
  temp_random_SILAC_2 = random_SILAC_Gene_2[,i]
  random_SILAC_cor[i] = cor(temp_random_SILAC_1,temp_random_SILAC_2, use= "pairwise.complete.obs",  method= c("pearson"))
}
random_pairs_cor = data.frame(Random_Gene_1 = colnames(random_TPM_Gene_1), Random_Gene_2 = colnames(random_TPM_Gene_2),
                               mRNA_correlation = random_TPM_cor, Protein_correlation = random_SILAC_cor)  

sig_RNA = wilcox.test(TPM_cor, random_TPM_cor)$p.value
sig_protein = wilcox.test(SILAC_cor, random_SILAC_cor)$p.value
sig_RNA_protein = wilcox.test(TPM_cor, SILAC_cor)$p.value

# Design the plot labels
mRNA_label = paste("mRNA m =", round(median(TPM_cor,na.rm=T),2))
mRNA_shuffled_label = paste("shuffled m =", round(median(random_TPM_cor,na.rm=T),2))
protein_label = paste("protein m =", round(median(SILAC_cor,na.rm=T),2))
protein_shuffled_label = paste("shuffled m =", round(median(random_SILAC_cor,na.rm=T),2))

sig_RNA_label = paste("P =", signif(sig_RNA, digits=1))
sig_protein_label = paste("P =", signif(sig_protein, digits=1))
sig_RNA_protein_label = paste("P =", signif(sig_RNA_protein, digits=1))

# Calculate the maximum bar height (as a help to set ylim) 
D1 = hist(gene_pairs_cor$Protein_correlation, breaks=seq(-1,1,0.1))$counts
D2 = hist(gene_pairs_cor$mRNA_correlation, breaks=seq(-1,1,0.1))$counts
D3 = hist(random_pairs_cor$Protein_correlation, breaks=seq(-1,1,0.1))$counts
D4 = hist(random_pairs_cor$mRNA_correlation, breaks=seq(-1,1,0.1))$counts
max_plot_height = max(D1, D2, D3, D4) * 1.02

# Set label positions in y dimension, evenly spaced and relative to the plot height
label_pos1 = max_plot_height*1.0
label_pos2 = max_plot_height*0.9
label_pos3 = max_plot_height*0.8

# Plot consecutive reaction pairs on mRNA and protein level vs random
p = ggplot(gene_pairs_cor)+
  xlim(-1,1)+
  ylim(0,max_plot_height)+
  xlab("Gene co-regulation [PCC]")+
  ylab("Gene pair count")+
  theme(panel.grid=element_blank(), panel.background=element_blank(), line=element_line(size=0.25),
        axis.text=element_text(size=6, colour="black"), axis.title=element_text(size=7),
        axis.line.x=element_line(colour="black"),axis.line.y=element_line(colour="black"),
        plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) 

pRNA = p+geom_histogram(aes(x=mRNA_correlation), binwidth=0.1, center=0.05, fill="#00A0BB", alpha=0.7)+
  geom_histogram(data=random_pairs_cor, aes(x=mRNA_correlation), binwidth=0.1, center=0.05, fill=NA, colour="grey30", size=0.25)+
  annotate("text", x=-0.98, y=label_pos1, label=mRNA_label, hjust=0, colour="#00A0BB", size=2.5)+
  annotate("text", x=-0.98, y=label_pos2, label=mRNA_shuffled_label, hjust=0, colour="grey30", size=2.5)+
  annotate("text", x=-0.98, y=label_pos3, label=sig_RNA_label, hjust=0, fontface="italic", size=2.5)

pProtein = p+geom_histogram(aes(x=Protein_correlation), binwidth=0.1, center=0.05, fill="#EC008C", alpha=0.7)+
  geom_histogram(data=random_pairs_cor, aes(x=Protein_correlation), binwidth=0.1, center=0.05, fill=NA, colour="grey30", size=0.25)+
  annotate("text", x=-0.98, y=label_pos1, label=protein_label, hjust=0, colour="#EC008C", size=2.5)+
  annotate("text", x=-0.98, y=label_pos2, label=protein_shuffled_label, hjust=0, colour="grey30", size=2.5)+
  annotate("text", x=-0.98, y=label_pos3, label=sig_protein_label, hjust=0, fontface="italic", size=2.5)+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank())

pRNA_Protein = p+geom_histogram(aes(x=mRNA_correlation), binwidth=0.1, center=0.05, fill="#00A0BB", alpha=0.7)+
  geom_histogram(aes(x=Protein_correlation), binwidth=0.1, center=0.05, fill="#EC008C", alpha=0.7)+
  annotate("text", x=-0.98, y=label_pos1, label=mRNA_label, hjust=0, colour="#00A0BB", size=2.5)+
  annotate("text", x=-0.98, y=label_pos2, label=protein_label, hjust=0, colour="#EC008C", size=2.5)+
  annotate("text", x=-0.98, y=label_pos3, label=sig_RNA_protein_label, hjust=0, fontface="italic", size=2.5)+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank())

pRNA_Protein

############# finished Figure 1C consecutive reactions #############

############# Figure 1C direct complexes (this is re-used code from above but with different test set) #############
# Load data
df = read.csv( file.path("data","mouse_SILAC_TPMs_log2_final_min8_features.csv"), stringsAsFactors=FALSE)
df = df[ !is.na( df[, "Ensembl"] ) ,]
rownames(df) = df[, "Ensembl"]
TPM = df[, grep("TPM_", colnames(df))]
SILAC = df[, grep("SILAC_", colnames(df))]
tTPM = t(TPM)
TPM_tissue_medians = apply(tTPM, 1, median, na.rm=TRUE)
tTPM_mn = sweep(tTPM, 1, TPM_tissue_medians, FUN="-")
tSILAC = t(SILAC)
SILAC_tissue_medians = apply(tSILAC, 1, median, na.rm=TRUE)
tSILAC_mn = sweep(tSILAC, 1, SILAC_tissue_medians, FUN="-")

# Load test data
all_test_pairs = read.csv( file.path("data","Reactome_corrs_dir_complexes_mouse.csv"), stringsAsFactors = FALSE)
is.na(all_test_pairs$Gene_1) = !(all_test_pairs$Gene_1 %in% rownames(df)) # turn NA if gene is NOT in our data
is.na(all_test_pairs$Gene_2) = !(all_test_pairs$Gene_2 %in% rownames(df)) # turn NA if gene is NOT in our data
our_gene_pairs = all_test_pairs[complete.cases(all_test_pairs),]   # remove rows with NA

# Compile matrices with input data for Gene_1 and Gene_2, respectively
TPM_Gene_1 = tTPM_mn[, match(our_gene_pairs$Gene_1, colnames(tTPM_mn))]
TPM_Gene_2 = tTPM_mn[, match(our_gene_pairs$Gene_2, colnames(tTPM_mn))]
SILAC_Gene_1 = tSILAC_mn[, match(our_gene_pairs$Gene_1, colnames(tSILAC_mn))]
SILAC_Gene_2 = tSILAC_mn[, match(our_gene_pairs$Gene_2, colnames(tSILAC_mn))]
TPM_cor = numeric()
SILAC_cor = numeric()

for (i in 1:nrow(our_gene_pairs)) {
  temp_TPM_1 = TPM_Gene_1[,i]
  temp_TPM_2 = TPM_Gene_2[,i]
  TPM_cor[i] = cor(temp_TPM_1,temp_TPM_2, use= "pairwise.complete.obs",  method= c("pearson"))
  temp_SILAC_1 = SILAC_Gene_1[,i]
  temp_SILAC_2 = SILAC_Gene_2[,i]
  SILAC_cor[i] = cor(temp_SILAC_1,temp_SILAC_2, use= "pairwise.complete.obs",  method= c("pearson"))
}
gene_pairs_cor = data.frame(Gene_1 = colnames(TPM_Gene_1), Gene_2 = colnames(TPM_Gene_2),
                            mRNA_correlation = TPM_cor, Protein_correlation = SILAC_cor)  
gene_pairs_cor = merge(gene_pairs_cor, our_gene_pairs, sort=FALSE)

# Create control by shuffling gene pairs

# Conditional shuffling here means that the test genes are paired randomly, but under the following conditions:
# (a) genes are never randomly paired with themselves
# (b) all randomly generated gene pairs are unique

# IMPORTANT: Re-execute this section if you receive an error message. The error occurs only if no more random pairs
# can be created due to the above restrictions.
gene_pool = c(our_gene_pairs$Gene_1, our_gene_pairs$Gene_2)             # Create a pool with all test genes

random_pairs_1 = character()                                            # Create empty vectors to fill up in loop
random_pairs_2 = character() 

while ( length(gene_pool) > 0) {
  
  Random_Gene_1 = sample( gene_pool ,1)                                   # Randomly choose a new Gene_1 
  
  previous_match_A = random_pairs_2[ random_pairs_1 == Random_Gene_1 ]  # When new Gene_1 was chosen as Random_gene_1 in previous iteration, which genes were paired as Random_gene_2?
  previous_match_B = random_pairs_1[ random_pairs_2 == Random_Gene_1 ]  # When new Gene_1 was chosen as Random_gene_2 in previous iteration, which genes were paired as Random_gene_1?
  previous_matches = c(previous_match_A, previous_match_B)              # All genes that were randomly matched to the current gene 1 in previous iterations
  previous_matches_indices = which( gene_pool %in% previous_matches )   # Indices of all genes that were previously matched to current gene 1
  
  self_indices = which( gene_pool == Random_Gene_1 )                    # Indices of all occurrences of Random_Gene_1, to avoid pairing it up with itself
  
  exclude = c(previous_matches_indices, self_indices)                   # Genes that must not be drawn as random gene 2
  
  Random_Gene_2 = sample( gene_pool[-exclude] ,1)                       # Randomly choose a new Gene_2 (with exceptions)
  
  random_pairs_1 = c(random_pairs_1, Random_Gene_1)                     # Add the randomly chosen genes to the appropriate vectors
  random_pairs_2 = c(random_pairs_2, Random_Gene_2)
  
  gene_pool = gene_pool[ -match( c(Random_Gene_1, Random_Gene_2) , gene_pool) ]    # Remove randomly chosen gene pair from the pool
}

random_gene_pairs = data.frame(Random_Gene_1=random_pairs_1, Random_Gene_2=random_pairs_2, stringsAsFactors = FALSE)    # Create dataframe with random gene pairs

# Compile matrices with input data for Random_Gene_1 and Random_Gene_2, respectively
random_TPM_Gene_1 = tTPM_mn[, match(random_gene_pairs$Random_Gene_1, colnames(tTPM_mn))]
random_TPM_Gene_2 = tTPM_mn[, match(random_gene_pairs$Random_Gene_2, colnames(tTPM_mn))]
random_SILAC_Gene_1 = tSILAC_mn[, match(random_gene_pairs$Random_Gene_1, colnames(tSILAC_mn))]
random_SILAC_Gene_2 = tSILAC_mn[, match(random_gene_pairs$Random_Gene_2, colnames(tSILAC_mn))]
random_TPM_cor = numeric()
random_SILAC_cor = numeric()

for (i in 1:nrow(random_gene_pairs)) {
  temp_random_TPM_1 = random_TPM_Gene_1[,i]
  temp_random_TPM_2 = random_TPM_Gene_2[,i]
  random_TPM_cor[i] = cor(temp_random_TPM_1,temp_random_TPM_2, use= "pairwise.complete.obs",  method= c("pearson"))
  temp_random_SILAC_1 = random_SILAC_Gene_1[,i]
  temp_random_SILAC_2 = random_SILAC_Gene_2[,i]
  random_SILAC_cor[i] = cor(temp_random_SILAC_1,temp_random_SILAC_2, use= "pairwise.complete.obs",  method= c("pearson"))
}
random_pairs_cor = data.frame(Random_Gene_1 = colnames(random_TPM_Gene_1), Random_Gene_2 = colnames(random_TPM_Gene_2),
                              mRNA_correlation = random_TPM_cor, Protein_correlation = random_SILAC_cor)  

sig_RNA = wilcox.test(TPM_cor, random_TPM_cor)$p.value
sig_protein = wilcox.test(SILAC_cor, random_SILAC_cor)$p.value
sig_RNA_protein = wilcox.test(TPM_cor, SILAC_cor)$p.value

# Design the plot labels
mRNA_label = paste("mRNA m =", round(median(TPM_cor,na.rm=T),2))
mRNA_shuffled_label = paste("shuffled m =", round(median(random_TPM_cor,na.rm=T),2))
protein_label = paste("protein m =", round(median(SILAC_cor,na.rm=T),2))
protein_shuffled_label = paste("shuffled m =", round(median(random_SILAC_cor,na.rm=T),2))

sig_RNA_label = paste("P =", signif(sig_RNA, digits=1))
sig_protein_label = paste("P =", signif(sig_protein, digits=1))
sig_RNA_protein_label = paste("P =", signif(sig_RNA_protein, digits=1))

# Calculate the maximum bar height (as a help to set ylim) 
D1 = hist(gene_pairs_cor$Protein_correlation, breaks=seq(-1,1,0.1))$counts
D2 = hist(gene_pairs_cor$mRNA_correlation, breaks=seq(-1,1,0.1))$counts
D3 = hist(random_pairs_cor$Protein_correlation, breaks=seq(-1,1,0.1))$counts
D4 = hist(random_pairs_cor$mRNA_correlation, breaks=seq(-1,1,0.1))$counts
max_plot_height = max(D1, D2, D3, D4) * 1.02

# Set label positions in y dimension, evenly spaced and relative to the plot height
label_pos1 = max_plot_height*1.0
label_pos2 = max_plot_height*0.9
label_pos3 = max_plot_height*0.8

# Plot consecutive reaction pairs on mRNA and protein level vs random
p = ggplot(gene_pairs_cor)+
  xlim(-1,1)+
  ylim(0,max_plot_height)+
  xlab("Gene co-regulation [PCC]")+
  ylab("Gene pair count")+
  theme(panel.grid=element_blank(), panel.background=element_blank(), line=element_line(size=0.25),
        axis.text=element_text(size=6, colour="black"), axis.title=element_text(size=7),
        axis.line.x=element_line(colour="black"),axis.line.y=element_line(colour="black"),
        plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) 

pRNA = p+geom_histogram(aes(x=mRNA_correlation), binwidth=0.1, center=0.05, fill="#00A0BB", alpha=0.7)+
  geom_histogram(data=random_pairs_cor, aes(x=mRNA_correlation), binwidth=0.1, center=0.05, fill=NA, colour="grey30", size=0.25)+
  annotate("text", x=-0.98, y=label_pos1, label=mRNA_label, hjust=0, colour="#00A0BB", size=2.5)+
  annotate("text", x=-0.98, y=label_pos2, label=mRNA_shuffled_label, hjust=0, colour="grey30", size=2.5)+
  annotate("text", x=-0.98, y=label_pos3, label=sig_RNA_label, hjust=0, fontface="italic", size=2.5)

pProtein = p+geom_histogram(aes(x=Protein_correlation), binwidth=0.1, center=0.05, fill="#EC008C", alpha=0.7)+
  geom_histogram(data=random_pairs_cor, aes(x=Protein_correlation), binwidth=0.1, center=0.05, fill=NA, colour="grey30", size=0.25)+
  annotate("text", x=-0.98, y=label_pos1, label=protein_label, hjust=0, colour="#EC008C", size=2.5)+
  annotate("text", x=-0.98, y=label_pos2, label=protein_shuffled_label, hjust=0, colour="grey30", size=2.5)+
  annotate("text", x=-0.98, y=label_pos3, label=sig_protein_label, hjust=0, fontface="italic", size=2.5)+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank())

pRNA_Protein = p+geom_histogram(aes(x=mRNA_correlation), binwidth=0.1, center=0.05, fill="#00A0BB", alpha=0.7)+
  geom_histogram(aes(x=Protein_correlation), binwidth=0.1, center=0.05, fill="#EC008C", alpha=0.7)+
  annotate("text", x=-0.98, y=label_pos1, label=mRNA_label, hjust=0, colour="#00A0BB", size=2.5)+
  annotate("text", x=-0.98, y=label_pos2, label=protein_label, hjust=0, colour="#EC008C", size=2.5)+
  annotate("text", x=-0.98, y=label_pos3, label=sig_RNA_protein_label, hjust=0, fontface="italic", size=2.5)+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank())

pRNA_Protein

############# finished Figure 1C direct complexes #############



############# Figure 1D #############
## Load mouse gene positions
genepos = fread( file.path("data","Mouse_gene_starts.csv") )

## Load and filter mRNA and protein level correlations
corrs = fread( file.path("data","Mouse_TPM_SILAC_correlations_pvals_sorted.csv"))
corrs = corrs[num_comm_obs_mRNA >= 12 & num_comm_obs_prot >= 12,] ## set minimum number of common observations
corrs$Gene_1_start = genepos$Gene_start[match(corrs$Gene_1,genepos$`Gene stable ID`)]
corrs$Gene_2_start = genepos$Gene_start[match(corrs$Gene_2,genepos$`Gene stable ID`)]
corrs$Gene_1_chr = genepos$`Chromosome/scaffold name`[match(corrs$Gene_1,genepos$`Gene stable ID`)]
corrs$Gene_2_chr = genepos$`Chromosome/scaffold name`[match(corrs$Gene_2,genepos$`Gene stable ID`)]
corrs = corrs[Gene_1_chr == Gene_2_chr,]
corrs = corrs[Gene_1 != Gene_2,]
corrs[,c("distance") := abs(Gene_1_start - Gene_2_start)]

## Divide genes into two groups, nearby genes (<= 50.000bp) and "others" (>50.000bp)
fifty = corrs[distance <= 5e4,]
others = corrs[distance > 5e4,]
fifty = fifty[,.(mRNA_correlation,prot_correlation,mRNA_BH_adj_P,prot_BH_adj_P)]
fifty$group = "close"
others = others[,.(mRNA_correlation,prot_correlation,mRNA_BH_adj_P,prot_BH_adj_P)]
others$group = "distant"
## Merge both groups for plotting and melt for long format
both= rbind(fifty,others)
both[,c("mRNA_BH_adj_P", "prot_BH_adj_P") := NULL]
both = melt(both,id.var="group")

## Calculate percentages
fifty_perc = data.table(c(nrow(fifty[mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05,])/nrow(fifty),
                          nrow(fifty[mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05 & prot_correlation > 0.5 & prot_BH_adj_P < 0.05,])/nrow(fifty)))
fifty_perc$group1 = c("mRNA","prot")
fifty_perc$group2 = "fifty"
other_perc = data.table(c(nrow(others[mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05,])/nrow(others),
                          nrow(others[mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05 & prot_correlation > 0.5 & prot_BH_adj_P < 0.05,])/nrow(others)))
other_perc$group1 = c("mRNA","prot")
other_perc$group2 = "others"
perc_both = rbind(fifty_perc,other_perc)

## Get counts instead of percentages and calculate significance using Fisher's exact test
close_mrna = nrow(fifty[mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05,])
close_prot = nrow(fifty[mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05 & prot_correlation > 0.5 & prot_BH_adj_P < 0.05,])
close_count = nrow(fifty)
distant_mrna = nrow(others[mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05,])
distant_prot = nrow(others[mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05 & prot_correlation > 0.5 & prot_BH_adj_P < 0.05,])
distant_count = nrow(others)

## Create 2 contigency tables for testing mRNA close vs distant, protein close vs distant
mrna_cont = matrix( c(close_mrna,close_count-close_mrna,distant_mrna,distant_count-distant_mrna), nrow=2, byrow=T)
prot_cont = matrix( c(close_prot,close_count-close_prot,distant_prot,distant_count-distant_prot), nrow=2, byrow=T)

## Test for significance
mrna_close_vs_dist_pval = chisq.test(mrna_cont)$p.value
prot_close_vs_dist_pval = chisq.test(prot_cont)$p.value

## Plot
ggplot(perc_both) + geom_bar(aes(x=group1,y=V1,fill=group2),stat="identity",position="dodge") + theme_bw() +
  ylab("Percentage of gene pairs") + xlab("") +
  scale_fill_discrete(labels=c(paste0("50kb_",nrow(fifty)),paste0("rest_",nrow(others)))) +
  ylim(0,0.15)

############# finished Figure 1D #############

############# Figure 1E #############
## Plot Chr17 and Chr2 correlation matricies
corrs = fread( file.path("data","Mouse_TPM_SILAC_correlations_pvals_sorted.csv") )
corrs = corrs[corrs$num_comm_obs_mRNA > 8 & corrs$num_comm_obs_prot > 8,] ## use 8 as minimum overlapping points to calculate correlations
corrs[,c("num_comm_obs_mRNA","num_comm_obs_prot") := NULL]

genestarts = fread(file.path("data","Mouse_gene_starts.csv"))
genestarts = genestarts[genestarts$`Chromosome/scaffold name` %in% c(1:22,"X","Y","MT"),]
genestarts = genestarts[genestarts$`Gene stable ID`%in% corrs$Gene_1 | genestarts$`Gene stable ID` %in% corrs$Gene_2 ,]
genestarts$`Chromosome/scaffold name` = paste("chr",genestarts$`Chromosome/scaffold name`,sep="")
genestarts[,Genomic_pos := rank(Gene_start),by=`Chromosome/scaffold name`]
corrs$Gene_1_chr = genestarts$`Chromosome/scaffold name`[match(corrs$Gene_1,genestarts$`Gene stable ID`)]
corrs$Gene_2_chr = genestarts$`Chromosome/scaffold name`[match(corrs$Gene_2,genestarts$`Gene stable ID`)]

## leave only same chr (not looking at correlations between chromosomes)
corrs = corrs[Gene_1_chr == Gene_2_chr,]
corrs$Gene_2_chr = NULL
colnames(corrs)[colnames(corrs) == "Gene_1_chr"] = "Chr"

## order genes by genomic position
corrs$Gene_1_start = genestarts$Gene_start[match(corrs$Gene_1,genestarts$`Gene stable ID`)]
corrs$Gene_2_start = genestarts$Gene_start[match(corrs$Gene_2,genestarts$`Gene stable ID`)]
corrs$Gene_1_pos = genestarts$Genomic_pos[match(corrs$Gene_1,genestarts$`Gene stable ID`)]
corrs$Gene_2_pos = genestarts$Genomic_pos[match(corrs$Gene_2,genestarts$`Gene stable ID`)]

## add second triangle 
corrs2 = copy(corrs)
corrs2[, c("Gene_2_pos", "Gene_1_pos") := .(Gene_1_pos, Gene_2_pos)] 
corrs2[, c("Gene_2", "Gene_1") := .(Gene_1, Gene_2)] 
corrs2[, c("Gene_2_start", "Gene_1_start") := .(Gene_1_start, Gene_2_start)] 
for_plotting = rbind(corrs,corrs2) ## now gene names are disconnected from their positions! NOT USE FOR OTHER ANALYSES
plotMatrix = function(selchr){
  selected = for_plotting[Chr==selchr,]
  selected = selected[order(selected$Gene_1_pos,selected$Gene_2_pos),]
  
  
  mplot = ggplotGrob(ggplot(selected,aes(Gene_1_pos,Gene_2_pos,fill=mRNA_correlation)) + geom_raster() + 
                       scale_fill_gradient2(low="#0080bb",mid="#eeeeee",high="#ff2500",midpoint = 0) + theme(panel.grid=element_blank(), panel.background=element_blank(), line=element_line(size=0.75),
                                                                                                             axis.text=element_text(size=10, colour="black"), axis.title=element_text(size=12),
                                                                                                             axis.line.x=element_line(colour="black"),axis.line.y=element_line(colour="black"),
                                                                                                             plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) + 
                       ggtitle(paste("mRNA correlation matrix for chromosome ",selchr,sep="")))
  
  pplot = ggplotGrob(ggplot(selected,aes(Gene_1_pos,Gene_2_pos,fill=prot_correlation)) + geom_raster() + 
                       scale_fill_gradient2(low="#0080bb",mid="#eeeeee",high="#ff2500",midpoint = 0) + theme(panel.grid=element_blank(), panel.background=element_blank(), line=element_line(size=0.75),
                                                                                                             axis.text=element_text(size=10, colour="black"), axis.title=element_text(size=12),
                                                                                                             axis.line.x=element_line(colour="black"),axis.line.y=element_line(colour="black"),
                                                                                                             plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) + 
                       ggtitle(paste("protein correlation matrix for chromosome ",selchr,sep="")))
  
  plot(mplot)
  plot(pplot)
}
plotMatrix("chr17")
plotMatrix("chr2")

############# finished Figure 1E #############

############# Figure 1F #############
## Plot Chr17 and Chr2 correlation profiles
genepos = fread( file.path("data","Mouse_gene_starts.csv"))
corrs = fread( file.path("data", "Mouse_TPM_SILAC_correlations_pvals_sorted.csv"))
corrs = corrs[num_comm_obs_mRNA >= 12 & num_comm_obs_prot >= 12,]
corrs$Gene_1_start = genepos$Gene_start[match(corrs$Gene_1,genepos$`Gene stable ID`)]
corrs$Gene_2_start = genepos$Gene_start[match(corrs$Gene_2,genepos$`Gene stable ID`)]
corrs$Gene_1_chr = genepos$`Chromosome/scaffold name`[match(corrs$Gene_1,genepos$`Gene stable ID`)]
corrs$Gene_2_chr = genepos$`Chromosome/scaffold name`[match(corrs$Gene_2,genepos$`Gene stable ID`)]
corrs = corrs[Gene_1_chr == Gene_2_chr,]
corrs = corrs[Gene_1 != Gene_2,]
corrs[,c("distance") := abs(Gene_1_start - Gene_2_start)]

chr17 = corrs[Gene_1_chr == "17",]
chr17 = chr17[,.(mRNA_correlation,prot_correlation,distance)]
chr17m = melt(chr17,id.vars = "distance")

chr17_plot = ggplot(chr17m) + 
  stat_smooth(aes(x=distance,y=value,color=variable),method="auto") + theme_bw() + 
  geom_hline(yintercept = 0, linetype="dashed",alpha=0.8) +
  theme(panel.grid=element_blank(), panel.background=element_blank(), line=element_line(size=0.25),
        axis.text=element_text(size=6, colour="black"), axis.title=element_text(size=7),
        axis.line = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.25),
        legend.position = "none") +
  coord_cartesian(xlim=c(0,6e7),ylim=c(-0.2,0.3)) 
chr17_plot

chr2 = corrs[Gene_1_chr == "2",]
chr2 = chr2[,.(mRNA_correlation,prot_correlation,distance)]
chr2m = melt(chr2,id.vars = "distance")

chr2m_plot = ggplot(chr2m) + 
  stat_smooth(aes(x=distance,y=value,color=variable),method="auto") + theme_bw() + 
  geom_hline(yintercept = 0, linetype="dashed",alpha=0.8) +
  theme(panel.grid=element_blank(), panel.background=element_blank(), line=element_line(size=0.25),
        axis.text=element_text(size=6, colour="black"), axis.title=element_text(size=7),
        axis.line = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.25),
        legend.position = "none") +
  coord_cartesian(xlim=c(0,6e7),ylim=c(-0.2,0.3)) 

chr2m_plot
############# finished Figure 1F #############