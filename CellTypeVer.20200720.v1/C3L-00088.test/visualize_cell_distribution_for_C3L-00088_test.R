#!/usr/bin/env Rscript --vanilla

library(optparse)
library(monocle)
library(ggplot2)

library(stringr)
library(dplyr)



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

cat("log file is at: ",paste0(out_path,"monocle_visualization.txt"),"\n")
sink(paste0(out_path,"monocle_visualization.txt"),append=TRUE)

obj_path <- opt$monocle_obj
if (!file.exists(obj_path)){
	cat(obj_path," does not exist!")
	sink()
	stop()
}


cat("reading monocle object...\n")
combined_subset_pseudotime <- readRDS(obj_path)

cat("plotting...\n")
pdf(paste0(out_path,"plot_cell_distribution_along_trajectory.pdf"),width=12,height=10)
        plot_cell_trajectory(combined_subset_pseudotime, color_by = "Naming",cell_size=0.3,cell_link_size=0.3,state_number_size=1)+facet_wrap(~State, nrow = 3)+theme(aspect.ratio=1)
        plot_cell_trajectory(combined_subset_pseudotime, color_by = "is_malignant",cell_size=0.3)+facet_wrap(~State, nrow = 1)+theme(aspect.ratio=1)
        plot_cell_trajectory(combined_subset_pseudotime, color_by = "State",cell_size=0.3)+theme(aspect.ratio=1)
dev.off()

sink()


