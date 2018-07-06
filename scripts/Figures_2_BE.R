## Script for generating Figures 2 B and E in the manuscript
## Authors: Piotr Grabowski and Georg Kustatscher, 05.2018
library(data.table)
library(reshape2)
library(stringr)
library(GenomicRanges)
library(gplots)
library(ggplot2)
library(plyr)

setwd("") ## set the working directory, the data will be loaded using relative paths

############# Figure 2B kmeans clustering of correlations #############
## Load expression data and correlations
data = fread( file.path(".","data","mouse_SILAC_TPMs_log2_final_min8_features.csv"))
corrs = fread( file.path(".","data","Mouse_TPM_SILAC_correlations_pvals_sorted.csv"))
corrs = corrs[corrs$num_comm_obs_mRNA > 8 & corrs$num_comm_obs_prot > 8,]
corrs[,c("num_comm_obs_mRNA","num_comm_obs_prot") := NULL]
## add the lower triangle of the matrix
corrs_second_h = copy(corrs)
corrs_second_h[, c("Gene_2", "Gene_1") := .(Gene_1, Gene_2)] 
corrs_both = rbind(corrs,corrs_second_h)
mrna_matrix = acast(corrs_both, formula = Gene_1 ~ Gene_2, value.var = "mRNA_correlation")
mrna_matrix = round(mrna_matrix,3)
diag(mrna_matrix) = 1
nans_cols = apply(mrna_matrix,2,function(x) sum(is.na(x)))
mrna_matrix = mrna_matrix[,nans_cols == 0]
## cluster data - check optimal k value
mrna_perc_explained = c()
for(k in c(2:10)){
  print(k)
  clustered = kmeans(mrna_matrix,k,nstart=3,iter.max = 20)
  perc_explained = clustered$betweenss/clustered$totss
  mrna_perc_explained = c(mrna_perc_explained,perc_explained)
}
qplot(c(2:10),mrna_perc_explained*100,geom="line",xlab="k value",main = "% variance explained in mRNA")

## repeat for protein
prot_matrix = acast(corrs_both, formula = Gene_1 ~ Gene_2, value.var = "prot_correlation")
prot_matrix = round(prot_matrix,3)
diag(prot_matrix) = 1

## NAs present in the matrix, remove features which have the most
nans_cols_p = apply(prot_matrix,2,function(x) sum(is.na(x)))
prot_matrix = prot_matrix[,nans_cols_p == 0]

## cluster data - check k value
prot_perc_explained = c()
for(k in c(2:10)){
  print(k)
  clustered = kmeans(prot_matrix,k,nstart=3,iter.max = 20)
  perc_explained = clustered$betweenss/clustered$totss
  prot_perc_explained = c(prot_perc_explained,perc_explained)
}

qplot(c(2:10),prot_perc_explained*100,geom="line",xlab="k value",main = "% variance explained in protein")

## optimal k for mRNA = 3, optimal k for protein = 5
clustered_rna = kmeans(mrna_matrix,3,nstart=3,iter.max = 20)
clustered_prot = kmeans(prot_matrix,5,nstart=3,iter.max = 20)
clustered_dt = data.table("Gene"=names(clustered_rna$cluster),"mRNA_cluster"=clustered_rna$cluster)
clustered_dt$protein_cluster = clustered_prot$cluster[match(clustered_dt$Gene,names(clustered_prot$cluster))]

## clusters assigned, plot the cluster matrix as in Fig 2B
dat = clustered_dt
dat$mRNA_cluster = factor(dat$mRNA_cluster,levels=c(1,3,2))
dat = dat[order(dat$mRNA_cluster),]
dat$mRNA_cluster_order = 1:nrow(dat)
dat = dat[order(dat$protein_cluster),]
dat$protein_cluster_order = 1:nrow(dat)
## Add cluster order to the correlations datatable
corrs_both$Gene_1_mRNA_cl_order = dat$mRNA_cluster_order[match(corrs_both$Gene_1,dat$Gene)]
corrs_both$Gene_2_mRNA_cl_order = dat$mRNA_cluster_order[match(corrs_both$Gene_2,dat$Gene)]
corrs_both$Gene_1_prot_cl_order = dat$protein_cluster_order[match(corrs_both$Gene_1,dat$Gene)]
corrs_both$Gene_2_prot_cl_order = dat$protein_cluster_order[match(corrs_both$Gene_2,dat$Gene)]

## plot mRNA matrix
corrs_both = corrs_both[order(Gene_1_mRNA_cl_order, Gene_2_mRNA_cl_order),]
ggplot(corrs_both) + geom_raster(aes(x=Gene_1_mRNA_cl_order,y=Gene_2_mRNA_cl_order,fill=mRNA_correlation),
                                 stat="identity",position="dodge") +
  scale_fill_gradient2(low="#0080bb",mid="#eeeeee",high="#ff2500",midpoint = 0) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) + guides(fill = "none")

## plot protein matrix
corrs_both = corrs_both[order(Gene_1_prot_cl_order, Gene_2_prot_cl_order),]
ggplot(corrs_both) + geom_raster(aes(x=Gene_1_prot_cl_order,y=Gene_2_prot_cl_order,fill=prot_correlation),
                                 stat="identity",position="dodge") +
  scale_fill_gradient2(low="#0080bb",mid="#eeeeee",high="#ff2500",midpoint = 0)

############# finished Figure 2B #############

############# Figure 2E epigenetic enrichments in clusters #############
# load epigenetic signal data (output from multiBigWigSummary) and gene positions
genepos = fread( file.path(".","data","Mouse_Biomart_gene_pos_29062017_mm10.csv"))
chips = fread( file.path(".","data","chipseq_avg_over_gene_bodies.tab"))

## metadata.tsv holds info about the used bigwig files with epigenetic data
file_meta = fread( file.path(".","data","metadata.tsv"))
file_meta$`File accession` = paste(file_meta$`File accession`,".bigWig",paste="")
file_meta$`File accession` = gsub(" ","",file_meta$`File accession`)
file_meta$df_name = paste(file_meta$`Biosample term name`,gsub("-mouse","",file_meta$`Experiment target`),sep="_")
colnames(chips) = gsub("'","",colnames(chips))
colnames(chips)[4:63] = file_meta$df_name[match(colnames(chips)[4:63],file_meta$`File accession`)]
chips$start = chips$start +1 ## go back to 1-based coords from GRanges object
chips$pos = paste(chips$chr,chips$start,chips$end,sep="_")
genepos$`Chromosome/scaffold name` = paste('chr',genepos$`Chromosome/scaffold name`,sep="")
genepos$`Chromosome/scaffold name`[genepos$`Chromosome/scaffold name` == "chrMT"] = "chrM"
genepos$pos = paste(genepos$`Chromosome/scaffold name`,genepos$`Gene start (bp)`,genepos$`Gene end (bp)`,sep="_")
chips$gene = genepos$`Gene stable ID`[match(chips$pos,genepos$pos)]
chips$pos = NULL
## remove tissues not present in our analysis
chips = chips[,-grep("testis|placenta|cortical plate|bone marrow",colnames(chips)),with=F]
marks = unique(sapply(str_split(colnames(chips)[4:49],"_"),"[",2)) ## 7 marks available

## add kmeans coregulation clusters
clusters = fread( file.path(".","data","protein_mrna_coregulation_clusters.csv"))
chips$mRNA_cluster = clusters$mRNA_cluster[match(chips$gene,clusters$Gene)]
chips$protein_cluster = clusters$protein_cluster[match(chips$gene,clusters$Gene)]
chips = chips[,-c("chr","start","end","gene")]

## change into an array, first get means of clusters
mrna_chips = chips[!is.na(chips$mRNA_cluster),]
mrna_chips = mrna_chips[,lapply(.SD,mean),by=mRNA_cluster,.SDcols = colnames(mrna_chips)[1:46]]

## calculate column means
mrna_colmeans = colMeans(mrna_chips[,2:ncol(mrna_chips)])
mrna_chips[,2:ncol(mrna_chips)] = sweep(mrna_chips[,2:ncol(mrna_chips)], 2, mrna_colmeans, `/`)
mrna_chips[,2:ncol(mrna_chips)] = log2(mrna_chips[,2:ncol(mrna_chips)])

## re-order data for plotting mark tis tis tis... mark2 tis tis tis...
## note: here cluster 2 is our new cluster3 (mito), it has to be re-ordered
mrna_chips$mRNA_cluster[2] = 2
mrna_chips$mRNA_cluster[3] = 3
mrna_chips_melt = melt(mrna_chips,id.vars = "mRNA_cluster")
mrna_chips_melt$tissue = sapply(str_split(mrna_chips_melt$variable,"_"),"[",1)
mrna_chips_melt$mark = sapply(str_split(mrna_chips_melt$variable,"_"),"[",2)
mrna_chips_melt$mark = factor(mrna_chips_melt$mark,levels=unique(mrna_chips_melt$mark))
mrna_chips_melt$tissue = factor(mrna_chips_melt$tissue, levels=unique(mrna_chips_melt$tissue))
mrna_chips_melt$variable = paste(mrna_chips_melt$mark,mrna_chips_melt$tissue,sep="_")
mrna_chips_melt = mrna_chips_melt[order(mark,tissue),]

## repeat for protein
## change into an array, first get means of clusters
prot_chips = chips[!is.na(chips$protein_cluster),]
prot_chips = prot_chips[,lapply(.SD,mean),by=protein_cluster,.SDcols = colnames(prot_chips)[1:46]]

## calculate column means
prot_colmeans = colMeans(prot_chips[,2:ncol(prot_chips)])
prot_chips[,2:ncol(prot_chips)] = sweep(prot_chips[,2:ncol(prot_chips)], 2, prot_colmeans, `/`)
prot_chips[,2:ncol(prot_chips)] = log2(prot_chips[,2:ncol(prot_chips)])

## re-order data for plotting mark tis tis tis... mark2 tis tis tis...
## note: here cluster 2 is our new cluster3 (mito), it has to be re-ordered
prot_chips_melt = melt(prot_chips,id.vars = "protein_cluster")
prot_chips_melt$tissue = sapply(str_split(prot_chips_melt$variable,"_"),"[",1)
prot_chips_melt$mark = sapply(str_split(prot_chips_melt$variable,"_"),"[",2)
prot_chips_melt$mark = factor(prot_chips_melt$mark,levels=unique(prot_chips_melt$mark))
prot_chips_melt$tissue = factor(prot_chips_melt$tissue, levels=unique(prot_chips_melt$tissue))
prot_chips_melt$variable = paste(prot_chips_melt$mark,prot_chips_melt$tissue,sep="_")
prot_chips_melt = prot_chips_melt[order(mark,tissue),]

## prepare colors for the plot
myBreaks = seq(-2,2,0.1)
my_palette = colorRampPalette(c("orange","orange", "white","steelblue","steelblue"))(n = length(myBreaks)-1)

## plot
mplot = ggplot(mrna_chips_melt,aes(variable,mRNA_cluster)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2(low="#0080bb",mid="#ffffff",high="#ff2500",midpoint = 0, limits=c(-0.9,0.9)) +
  #theme_classic(base_size = 8) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  coord_fixed()

pplot = ggplot(prot_chips_melt,aes(variable,protein_cluster)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2(low="#0080bb",mid="#ffffff",high="#ff2500",midpoint = 0, limits=c(-0.9,0.9)) +
  #theme_classic(base_size = 8) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  coord_fixed()

mplot
pplot

## calculate coefficients of variation and plot the dots used in the figure
cv = function(vec){
  covar = abs(round(sd(vec)/mean(vec),1))
  covar
}
prot_cvs = prot_chips_melt[,cv(value),by=c("protein_cluster","mark")]
prot_cvs = prot_cvs[order(mark,protein_cluster),]
mrna_cvs = mrna_chips_melt[,cv(value),by=c("mRNA_cluster","mark")] ## the clusters 2/3 were already re-ordered
mrna_cvs$cluster = paste("T",mrna_cvs$mRNA_cluster,sep="")
prot_cvs$cluster = paste("P",prot_cvs$protein_cluster,sep="")

mrna_cvs$mRNA_cluster = NULL
prot_cvs$protein_cluster = NULL

both = rbind(mrna_cvs,prot_cvs)
both$V1 = both$V1 + 0.01

both$cluster = factor(both$cluster,levels=c("T1","T2","T3","P1","P2","P3","P4","P5"))
both$mark = factor(both$mark,levels=c("H3K27ac","H3K36me3","H3K4me1","H3K4me3","H3K9ac"))
both = both[!is.na(both$mark),]
both = both[order(cluster,mark),]
CV_dots = ggplot(both) + 
  geom_point(aes(x=mark,y=cluster,color=log2(V1)),size=5) + 
  scale_color_gradient2(low="white",high="#601E9E") +
  theme_classic(base_size = 6)
CV_dots

############# finished Figure 2E #############