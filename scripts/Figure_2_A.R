## A script to obtain plots about the 5 post-transcriptional mechanisms and arrange them aligned for 
## publication (Figure 2A)
## Author: Piotr Grabowski, 2018 Rappsilber lab
library(cowplot)

setwd("") ## set the working directory, the data will be loaded using relative paths

## Run the 4 scripts to get plot objects
source(file.path(".","scripts","miRNA_buffering.R"))
source(file.path(".","scripts","coding_length_buffering.R"))
source(file.path(".","scripts","ribo_profiling_buffering.R"))
source(file.path(".","scripts","NEDs_buffering.R"))
source(file.path(".","scripts","prot_translation_rates_buffering.R"))

## Clean workspace, keep the 4 final plots
rm(list=setdiff(ls(), c("cds_buff_plot","cod_len_plot","ribo_profile_plot","ned_plot", "sim_tl_rates_plot")))

## Remove legends from each
cod_len_plot = cod_len_plot + theme(legend.position="none", text = element_text(size=6), 
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank())
mirna_plot_cds = cds_buff_plot + theme(legend.position="none", text = element_text(size=6), 
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ribo_profile_plot = ribo_profile_plot + theme(legend.position="none", text = element_text(size=6),
                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ned_plot = ned_plot + theme(legend.position="none", text = element_text(size=6),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sim_tl_rates_plot = sim_tl_rates_plot + theme(legend.position="none", text = element_text(size=6),
                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Figure 2A
plot_grid(cod_len_plot,mirna_plot_cds,ribo_profile_plot,sim_tl_rates_plot,ned_plot,
                 align = "h", axis = b,
                 ncol=5, label_size = 8, label_fontfamily = "arial")
