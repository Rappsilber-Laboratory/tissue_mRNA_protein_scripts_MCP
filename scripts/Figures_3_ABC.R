## Script for generating Figures 3 A,B,C in the manuscript
## Authors: Piotr Grabowski and Georg Kustatscher, 05.2018

library(psych)
library(ggplot2)
library(data.table)

setwd("") ## set the working directory, the data will be loaded using relative paths

## sef-error proof corr.test:
corr_mouse_test = function(x, y = NULL, use = "pairwise", method = "pearson", 
                            adjust = "holm", alpha = 0.05) 
{
  cl = match.call()
  if (is.null(y)) {
    r = cor(x, use = use, method = method)
    sym = TRUE
    n = t(!is.na(x)) %*% (!is.na(x))
  }
  else {
    r = cor(x, y, use = use, method = method)
    sym = FALSE
    n = t(!is.na(x)) %*% (!is.na(y))
  }
  if ((use == "complete") | (min(n) == max(n))) 
    n = min(n)
  t = (r * sqrt(n - 2))/sqrt(1 - r^2)
  p = 2 * (1 - pt(abs(t), (n - 2)))
  se = sqrt((1 - r * r)/(n - 2))
  nvar = ncol(r)
  p[p > 1] = 1
  if (adjust != "none") {
    if (is.null(y)) {
      lp = upper.tri(p)
      pa = p[lp]
      pa = p.adjust(pa, adjust)
      p[upper.tri(p, diag = FALSE)] = pa
    }
    else {
      p[] = p.adjust(p, adjust)
    }
  }
  z = fisherz(r[lower.tri(r)])
  
  ci = NULL
  result = list(r = r, n = n, t = t, p = p, se = se, 
                 adjust = adjust, sym = sym, 
                 Call = cl)
  class(result) = c("psych", "corr.test")
  return(result)
}

############# Figure 3A expression variability on mRNA and protein levels #############
# Load and prep the input data
df = read.csv( file.path(".","data","mouse_SILAC_TPMs_log2_final_min8_features.csv"), stringsAsFactors = FALSE)
rownames(df) = df[, "Ensembl"]
RPKM = df[, grep("TPM_", colnames(df))] ## it's transcripts-per-million in the mouse data
SILAC = df[, grep("SILAC_", colnames(df))]
# Transpose the matrixes and sweep out row-medians to avoid artificial correlations
tRPKM = t(RPKM)
RPKM_people_medians = apply(tRPKM, 1, median, na.rm=TRUE)
tRPKM_mn = sweep(tRPKM, 1, RPKM_people_medians, FUN="-")

tSILAC = t(SILAC)
SILAC_people_medians = apply(tSILAC, 1, median, na.rm=TRUE)
tSILAC_mn = sweep(tSILAC, 1, SILAC_people_medians, FUN="-")

# Get all possible mRNA correlations
tRPKM_cor = corr_mouse_test(tRPKM_mn, use = "pairwise", method = c("pearson"), adjust="BH")    # Get correlation matrix with P values
mRNA_correlations = tRPKM_cor$r
mRNA_pvalues      = tRPKM_cor$p
mRNA_correlations[ lower.tri(mRNA_correlations, diag=TRUE) ] = NA     # Set lower tri and diag to NA
mRNA_pvalues[      lower.tri(mRNA_pvalues,      diag=TRUE) ] = NA     # Set lower tri and diag to NA (contains non-adjusted p-values)
mRNA_correlations = melt(mRNA_correlations)
mRNA_pvalues      = melt(mRNA_pvalues)
colnames(mRNA_correlations) = c("Gene_1", "Gene_2", "mRNA_correlation")
colnames(mRNA_pvalues)      = c("Gene_1", "Gene_2", "mRNA_BH_adj_P")
mRNA_correlations = as.data.table(mRNA_correlations)
mRNA_pvalues      = as.data.table(mRNA_pvalues)
setkey(mRNA_correlations, Gene_1, Gene_2)
setkey(mRNA_pvalues,      Gene_1, Gene_2)
mRNA_DT = merge(mRNA_correlations, mRNA_pvalues)
mRNA_DT = mRNA_DT[ complete.cases(mRNA_DT) ]

# Get all possible protein correlations
tSILAC_cor = corr_mouse_test(tSILAC_mn, use = "pairwise", method = c("pearson"), adjust="BH")    # Get correlation matrix with P values
prot_correlations = tSILAC_cor$r
prot_pvalues      = tSILAC_cor$p
prot_correlations[ lower.tri(prot_correlations, diag=TRUE) ] = NA     # Set lower tri and diag to NA
prot_pvalues[      lower.tri(prot_pvalues,      diag=TRUE) ] = NA     # Set lower tri and diag to NA (contains non-adjusted p-values)
prot_correlations = melt(prot_correlations)
prot_pvalues      = melt(prot_pvalues)
colnames(prot_correlations) = c("Gene_1", "Gene_2", "prot_correlation")
colnames(prot_pvalues)      = c("Gene_1", "Gene_2", "prot_BH_adj_P")
prot_correlations = as.data.table(prot_correlations)
prot_pvalues      = as.data.table(prot_pvalues)
setkey(prot_correlations, Gene_1, Gene_2)
setkey(prot_pvalues,      Gene_1, Gene_2)
prot_DT = merge(prot_correlations, prot_pvalues)
prot_DT = prot_DT[ complete.cases(prot_DT) ]

# Merge the mRNA and protein tables
DT = merge(mRNA_DT, prot_DT)

# Load genome positions
Genes = fread( file.path(".","data","Mouse_Biomart_gene_pos_29062017_mm10.csv") )
colnames(Genes) = gsub(" ", ".", colnames(Genes), fixed = TRUE)
colnames(Genes) = gsub("(", ".", colnames(Genes), fixed = TRUE)
colnames(Genes) = gsub(")", ".", colnames(Genes), fixed = TRUE)
Genes$Gene_TSS = ifelse(Genes$Strand == 1, Genes$Gene.start..bp., Genes$Gene.end..bp.)   # Determine Gene TSSs
Genes$Gene_end = ifelse(Genes$Strand == 1, Genes$Gene.end..bp.,  Genes$Gene.start..bp.)  # Determine end of gene
setkey(Genes, Gene.stable.ID, Gene_TSS)

# For gene pairs from the same chromosome, calculate distance between TSSs
DT = DT[ (Gene_1 %in% Genes$Gene.stable.ID) & (Gene_2 %in% Genes$Gene.stable.ID) ]   # Keep only pairs for which we do have gene positions
DT = merge(DT, Genes, by.x="Gene_1", by.y="Gene.stable.ID", all.x = TRUE)             # Add genomic annotation for Gene 1
DT = merge(DT, Genes, by.x="Gene_2", by.y="Gene.stable.ID", all.x = TRUE)             # Add genomic annodation for Gene 2
DT = DT[DT$`Chromosome/scaffold.name.x` == DT$`Chromosome/scaffold.name.y` ,]
DT[, Distance := abs(DT$Gene_TSS.x-DT$Gene_TSS.y) ]   # For pairs on same chromosome, calculate distance between their TSSs
DT = DT[, .(Gene_1, Gene_2, mRNA_correlation, mRNA_BH_adj_P, prot_correlation, prot_BH_adj_P, Distance)]  # Keep relevant columns only

## Calculate gene expression noise (coefficient of variation)
# Re-load and prep the input data
df = read.csv( file.path(".",,"data","mouse_SILAC_TPMs_log2_final_min8_features.csv"), stringsAsFactors = FALSE)
rownames(df) = df[, "Ensembl"]
RPKM = df[, grep("TPM_", colnames(df))] ## it's transcripts-per-million in the mouse data
SILAC = df[, grep("SILAC_", colnames(df))]

tRPKM = t(RPKM)
RPKM_people_medians = apply(tRPKM, 1, median, na.rm=TRUE)
tRPKM_mn = sweep(tRPKM, 1, RPKM_people_medians, FUN="-")

tSILAC = t(SILAC)
SILAC_people_medians = apply(tSILAC, 1, median, na.rm=TRUE)
tSILAC_mn = sweep(tSILAC, 1, SILAC_people_medians, FUN="-")

# Shift ratios to avoid division by zero
tRPKM_mn  = tRPKM_mn+10      # This does not affect SD but shifts all means away from zero into positive values
tSILAC_mn = tSILAC_mn+10     # (avoiding negative noise as well as inflating noise by division with mean values close to zero)

# Calculate mRNA noise
mRNA_sd = apply(tRPKM_mn, 2, sd  , na.rm = TRUE)
mRNA_me = apply(tRPKM_mn, 2, mean, na.rm = TRUE)
mRNA_CV = mRNA_sd / mRNA_me

# Calculate protein noise
prot_sd = apply(tSILAC_mn, 2, sd  , na.rm = TRUE)
prot_me = apply(tSILAC_mn, 2, mean, na.rm = TRUE)
prot_CV = prot_sd / prot_me

# Combine data
noise = merge( as.data.frame(mRNA_CV) , as.data.frame(prot_CV) , by="row.names")
colnames(noise)[1] = "Gene_ID"
noise = as.data.table(noise)

#### Calculate local gene density for all human genes - linear contacts (bp) ####
intra = DT[ !is.na(Distance) ]                                       # Get gene pairs that are on the same chromosome
intra_genes = unique( intra[, c(Gene_1, Gene_2)] )                   # And the names of all the genes
intra_mean_dist_to5NN = data.table()                                 # Initialise result table
pb = txtProgressBar(min = 0, max = length(intra_genes), style = 3)   # Initialise progress bar

for(i in intra_genes){
  temp_dist = intra[ Gene_1 == i | Gene_2 == i, Distance]     # Get the distances (bp) to all other genes on chromosome
  temp_mean = mean( head( sort(temp_dist), n = 3))            # Get the mean distance of the 3 nearest genes
  intra_mean_dist_to5NN = rbind( intra_mean_dist_to5NN,       # Combine gene IDs and mean distances into a table
                                  data.table(Gene_ID = i,
                                             intra_mean_dist_to5NN = temp_mean))
  setTxtProgressBar(pb, intra_mean_dist_to5NN[,.N])
}

rm( list = ls()[! ls() %in% c("DT", "noise", "intra_mean_dist_to5NN")] )   # Clear workspace

# Combine noise levels and gene densities 
setkey(noise                   , Gene_ID)
setkey(intra_mean_dist_to5NN   , Gene_ID)
noise = merge(noise, intra_mean_dist_to5NN)

# Noise reduction with clustering
noise[,    intra_rank := frank(intra_mean_dist_to5NN) ]     # Ranking allows me to subset the most and least clustered genes afterwards

noise[ intra_rank <= ceiling(nrow(noise) * 0.05)    , intra := "most_5pc" ]      # Define the top and least categories   
noise[ intra_rank >= ceiling(nrow(noise) * 0.95)   , intra := "least_5pc" ]        

# Significance of mRNA expression noise reduction
wilcox.test( noise[ intra == "most_5pc", mRNA_CV ] , noise[ intra == "least_5pc", mRNA_CV ] )$p.value

# Significance of protein expression noise reduction
wilcox.test( noise[ intra == "most_5pc", prot_CV ] , noise[ intra == "least_5pc", prot_CV ] )$p.value

# Combine data for plotting
boxplot_data = rbind( 
  noise[ !is.na(intra)    , .(Gene_ID, value = mRNA_CV, type = intra, annotation = "intra mRNA")  ]       ,
  noise[ !is.na(intra)    , .(Gene_ID, value = prot_CV, type = intra, annotation = "intra prot")  ]        )

boxplot_data$type = factor(boxplot_data$type, levels = c("most_5pc", "least_5pc"))
boxplot_data$annotation = factor(boxplot_data$annotation, levels = c("intra mRNA","intra prot"))

p2 = ggplot(boxplot_data, aes(x = annotation, y = value, fill = type))+
  geom_boxplot(outlier.colour = NA, notch = T, size=0.25)+
  #coord_cartesian(ylim = c(0.005,0.09))+
  scale_fill_manual(values = c("darkorange1", "deepskyblue3"))+
  ylab("Gene expression noise (CV)")+
  theme(panel.grid=element_blank(), panel.background=element_blank(), line=element_line(size=0.25),
        axis.text=element_text(size=6, colour="black"), axis.title.y=element_text(size=7),
        axis.line.x=element_line(colour="black"),axis.line.y=element_line(colour="black"),
        panel.border = element_rect(fill=NA, colour="black", size=0.25), axis.line = element_blank(),
        axis.title.x = element_blank(), legend.position = "top", legend.title = element_blank(),
        legend.box.background = element_blank(), legend.key = element_blank())

p2

############# finished Figure 3A #############

############# Figure 3B epigenetic similarity vs mRNA and protein co-expression #############
epi_dists = fread( file.path(".","data","Mahalanobis_dist.csv") )
epi_dists[, Gene_1_sorted := ifelse(Gene_1 < Gene_2, Gene_1, Gene_2) ]                                # Sort the IDs row-wise
epi_dists[, Gene_2_sorted := ifelse(Gene_1 > Gene_2, Gene_1, Gene_2) ]    
lin_dists_corrs = fread( file.path(".","data","Mouse_TPM_SILAC_correlations_pvals_sorted.csv"))
## get linear separation for same chromosome
same_chr_separation = lin_dists_corrs[Gene_1_chr == Gene_2_chr,]
same_chr_separation[,separation := abs(Gene_1_start - Gene_2_start)]
## add epigenetic similarity measure
same_chr_separation$mahal_dist = epi_dists$maha_dist[match(paste(same_chr_separation$Gene_1,same_chr_separation$Gene_2),
                                                           paste(epi_dists$Gene_1_sorted,epi_dists$Gene_2_sorted))]

same_chr_separation[,epi_similarity := 1/same_chr_separation$mahal_dist]
same_chr_separation = same_chr_separation[epi_similarity < 0.2,]

dist_pccs = ggplot(same_chr_separation, aes(x = epi_similarity))+         
  geom_smooth(aes(y = mRNA_correlation), colour="#00A0BB", alpha=0.7)+
  geom_smooth(aes(y = prot_correlation), colour="#EC008C", alpha=0.7)+
  coord_cartesian( ylim = c(-0.3, 0.6))+
  geom_hline(yintercept = 0, linetype="dashed", colour="grey60", size=0.25)+
  ylab("Gene co-regulation [PCC]")+
  xlab("Epigenetic similarity \n[1/Mahalanobis distance]")+
  theme(panel.grid=element_blank(), panel.background=element_blank(), line=element_line(size=0.25),
        axis.text=element_text(size=6, colour="black"), axis.title=element_text(size=7),
        axis.line = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.25)) 

dist_pccs

############# finished Figure 3B #############

############# Figure 3C percentages of coregulated pairs as a function of epigenetic similarity and linear distances #############
# Bin INTRA gene pairs into groups with (approximately) equal numbers of observations,
numbins = 10
same_chr_separation[, linear_bins := cut_number(1/separation, numbins, labels = FALSE)]
same_chr_separation[, epi_bins := cut_number(epi_similarity, numbins, labels = FALSE)]

## Label correlated pairs
same_chr_separation[ mRNA_correlation > 0.5 & mRNA_BH_adj_P < 0.05 , mRNA_coregulated := "yes" ]
same_chr_separation[ prot_correlation > 0.5 & prot_BH_adj_P < 0.05 , prot_coregulated := "yes" ]

# Get percentage of co-regulated mRNAs / proteins per combination (sector) and plot results
Pintra = same_chr_separation[, .(pc_coreg_mRNA = sum(mRNA_coregulated == "yes", na.rm = TRUE)/.N*100,
                                  pc_coreg_prot = sum(prot_coregulated == "yes", na.rm = TRUE)/.N*100,
                                  N_coreg_mRNA = sum(mRNA_coregulated == "yes", na.rm = TRUE),
                                  N_coreg_prot = sum(prot_coregulated == "yes", na.rm = TRUE),
                                  N_sector = .N) , 
                              by = .(linear_bins, epi_bins) ]

## Plot
pIntra_mRNA = ggplot(Pintra, aes(x = linear_bins, y = epi_bins, colour = pc_coreg_mRNA, size = N_sector))+
  geom_point()+
  scale_size(breaks=c(2000,2500,3000), range=c(4.8,7))+
  xlab("binned by decreasing linear distance")+
  ylab("binned by increasing epigenetic similarity")+
  scale_colour_gradient(low="darkblue", high="deeppink1", breaks = c(5,10,15,20,25,30))+
  annotate(geom="rect", xmin = c(0.5, 0.5, 9.5, 9.5), xmax = c(1.5, 1.5, 10.5, 10.5),
           ymin = c(0.5, 9.5, 0.5, 9.5), ymax = c(1.5, 10.5, 1.5, 10.5), fill=NA, colour="black", size=0.25)+
  theme(panel.background = element_blank(), panel.grid=element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=0.25),
        axis.ticks = element_blank(), axis.text = element_blank(), legend.title=element_text(size=7),
        legend.text=element_text(size=6), strip.background = element_rect(colour="black", size=0.25), strip.text = element_text(size=7),
        axis.title = element_text(size=7))
pIntra_mRNA


pIntra_prot = ggplot(Pintra, aes(x = linear_bins, y = epi_bins, colour = pc_coreg_prot))+
  geom_point(size=2.5)+
  xlab("binned by decreasing linear distance")+
  ylab("binned by increasing epigenetic similarity")+
  scale_colour_gradient(low="darkblue", high="deeppink1", breaks = c(4,8,12,16), limits = c(0.8, 17.5))+
  theme(panel.background = element_blank(), panel.grid=element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=0.25),
        axis.ticks = element_blank(), axis.text = element_blank(), legend.title=element_text(size=7),
        legend.text=element_text(size=6), strip.background = element_rect(colour="black", size=0.25), strip.text = element_text(size=7),
        axis.title = element_text(size=7),legend.position = "none")
pIntra_prot

############# finished Figure 3C #############