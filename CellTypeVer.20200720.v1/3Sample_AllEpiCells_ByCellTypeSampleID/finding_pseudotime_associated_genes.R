## !! DO NOT RUN THIS IF THE MONOCLE OBJECT IS NOT UPDATED; TAKES A HUGE AMOUNT OF MEMORY AND TIME


library(optparse)
library(monocle)
library(stringr)
library(dplyr)

option_list = list(
	make_option(c("--monocle_obj"),
		type="character",
		default=NA,
		help="path to the monocle object to be analyzed",
		metavar="characer"),
	make_option(c("-o","--out_path"),
		type="character",
		default="./",
		help="path to the output path",
		metavar="character"),
	make_option(c("-b","--branch_point"),
		type="character",
		default=NA,
		help="specify which branch point to use to analyze branch-determining genes; when multiple branch points are needed, use ',' to separate. For example '1,2,3' ",
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
dir.create(out_path,showWarnings=TRUE,recursive=TRUE)

cat("log file is at: ",paste0(out_path,"monocle_DE_analysis.txt"),"\n")
sink(paste0(out_path,"monocle_DE_analysis.txt"),append=TRUE)


obj_path <- opt$monocle_obj
if (!file.exists(obj_path)){
        cat(obj_path," does not exist!")
        sink()
        stop()
}

cat("reading monocle object...\n")
combined_subset_pseudotime <- readRDS(obj_path)

branch_points <- opt$branch_point %>% strsplit(",") %>% unlist %>% as.numeric
cat("branch points to be analyzed:")
branch_points
cat("\n")

## find genes that change along pseudotime
cat("finding DE genes across pseudotime...\n")
diff_test_res <- differentialGeneTest(combined_subset_pseudotime,fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.table(diff_test_res,paste0(out_path,"DE_genes_ModelBy_Pseudotime.txt"),sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)

## find genes that determine the branch
if (!is.na(branch_points)){
	BEAM_res_list <- list()
	for (branch_point in branch_points){
		cat("finding branch point determining genes for branch point ",branch_point,"\n")
		BEAM_res_list[[branch_point]] <- BEAM(combined_subset_pseudotime, branch_point = branch_point, cores = 10)
		BEAM_res_list[[branch_point]] <- BEAM_res_list[[branch_point]][order(BEAM_res_list[[branch_point]]$qval),]
		BEAM_res_list[[branch_point]] <- BEAM_res_list[[branch_point]][,c("gene_short_name", "pval", "qval")]
		BEAM_res_list[[branch_point]]$branch_point <- branch_point
	}
	BEAM_res <- do.call("rbind",BEAM_res_list)
	write.table(BEAM_res,paste0(out_path,"branch_determining_genes.txt"),sep="\t",row.names=FALSE,quote=FALSE)
}

cat("done")
sink()
