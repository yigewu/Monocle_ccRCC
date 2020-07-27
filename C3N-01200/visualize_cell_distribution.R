#!/usr/bin/env Rscript --vanilla

# Set Run-Specific Parameters ---------------------------------------------

# Set Fixed Parameters -------------------------------------------------------------------
# load required libraries
library(optparse)
library(monocle)
library(ggplot2)
library(stringr)
library(dplyr)

# create user options
option_list = list(
	make_option(c("-o","--out_path"),
		type="character",
		default="./",
		help="path to the output path",
		metavar="character"),
	make_option(c("--monocle_obj"),
		type="character",
		default=NA,
		help="path to the monocle object to be visualized",
		metavar="character")
);

# read in initial arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## assign arguments to values
out_path <- opt$out_path
if (str_sub(out_path,"-1")!="/"){
  out_path<-paste0(out_path,"/")
}

## write log file
cat("log file is at: ",paste0(out_path,"monocle_visualization.txt"),"\n")
sink(paste0(out_path,"monocle_visualization.txt"),append=TRUE)

# Input Monocle Object if It Already Exists--------------------------------------------------
obj_path <- opt$monocle_obj
if (!file.exists(obj_path)){
	cat(obj_path," does not exist!")
	sink()
	stop()
}

cat("reading monocle object...\n")
combined_subset_pseudotime <- readRDS(obj_path)

# Plottting--------------------------------------------------
cat("plotting...\n")
pdf(paste0(out_path,"cell_distribution_along_trajectory.pdf"),width=12,height=10)
        plot_cell_trajectory(combined_subset_pseudotime, color_by = "Aliquot.snRNA.WU",cell_size=0.3,cell_link_size=0.3,state_number_size=1)+facet_wrap(~State, nrow = 3)+theme(aspect.ratio=1)
        plot_cell_trajectory(combined_subset_pseudotime, color_by = "Cell_type.detailed",cell_size=0.3)+facet_wrap(~State, nrow = 1)+theme(aspect.ratio=1)
        plot_cell_trajectory(combined_subset_pseudotime, color_by = "State",cell_size=0.3)+theme(aspect.ratio=1)
dev.off()
cat("Finished Plotting!\n")
sink()


