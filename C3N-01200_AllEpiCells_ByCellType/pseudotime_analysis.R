#!/usr/bin/env Rscript --vanilla


# Set Run-Specific Parameters ---------------------------------------------
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
path_srat <- "./Resources/snRNA_Processed_Data/Merged_Seurat_Objects/C3N-01200.Tumor_Segments.Merged.20200319.v1.RDS"
path_idmetadata <- "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv"
# path_bc2celltype <- "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells_with_patch/20200720.v1/31AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv"
path_bc2celltype <- "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200811.v1/31Aliquot.Barcode2CellType.20200811.v1.tsv"

# Set Fixed Parameters -------------------------------------------------------------------
# load required libraries
library(optparse) ## installed
library(Seurat) ## installed
library(monocle) ## installed
library(tibble) ## installed
library(dplyr) ## installed
library(stringr) ## installed
library(grid)

## create function
order_cells <- function(monocle_obj,ordering_genes,visualization_output,reverse=FALSE){
  monocle_obj <- setOrderingFilter(monocle_obj, ordering_genes)
  pdf(visualization_output)
  p <- plot_ordering_genes(monocle_obj)
  print(p)
  dev.off()
  
  ### reduce data dimensionality
  monocle_obj <- reduceDimension(monocle_obj, max_components = 2,method = 'DDRTree')
  
  ### order cells along the trajectory
  monocle_obj <- orderCells(monocle_obj,reverse=reverse)
  
  return (monocle_obj)
}

# create user options
option_list = list(
  make_option(c("--FormulaStr"),
              type="character",
              default=NA,
              help="specify what features are used for testing differential expression (the result of which would be used for ordering)",
              metavar="character"),
  make_option(c("-q", "--q_val"),
              type="character",
              default="0.01",
              help="q value cutoff for selecting features used for ordeinrg",
              metavar="character"),
  make_option(c("-o","--out_path"),
              type="character",
              default="./",
              help="path to the output path",
              metavar="character"),
  make_option(c("--force.reprocess"),
              type="logical",
              default=FALSE,
              help="whether to reprocess when the monocle object already exists",
              metavar="character"),
  make_option(c("--reverse"),
              type="logical",
              default=FALSE,
              help="whether to reverse the trajectory when the monocle object",
              metavar="character"),
  make_option(c("--root_state"),
              type="integer",
              default=NA,
              help="specify which state to set as the root state. Only valid when the monocle object already exists",
              metavar="integer")
);


# read in initial arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## assign arguments to values
out_path <- opt$out_path
if (str_sub(out_path,"-1")!="/"){
  out_path<-paste0(out_path,"/")
}

formula_str <- opt$FormulaStr

q_val_cutoff <- opt$q_val %>% as.numeric


## write log file
dir.create(out_path,recursive=TRUE,showWarnings = TRUE)
cat(paste0("Processing log is at \n", out_path, "monocle_processing.txt\n"))
sink(paste0(out_path, "monocle_processing.txt"),append=TRUE)
cat ("differential expression is based on: ",formula_str,"\n")
cat ("q value cutoff is:",q_val_cutoff %>% as.character,"\n")
cat ("output path is:",out_path,"\n")


# Input Monocle Object if It Already Exists, re-call OrderCells to define the "exact" root state --------------------------------------------------
## whether the monocle object already exists and whether it needs to obe re-processed
obj_output_path <- paste0(out_path,"combined_subset_pseudotime_qval_",q_val_cutoff,".rds")
if (file.exists(obj_output_path) & opt$force.reprocess==FALSE){
  cat (obj_output_path," already exists. Will not do the reprocess!\n")
  combined_subset_pseudotime <- readRDS(obj_output_path)
  root_state <- opt$root_state
  if (is.na(root_state)){
    root_state<-1
    cat("root state is not specified! Set 1 instead...\n")
  }
  reverse <- opt$reverse
  cat ("root  state is set as ",root_state,"\n")
  cat ("reverse trajectory: ",reverse,"\n")
  
  ## re-call OrderCells to define the "exact" root state
  combined_subset_pseudotime <- orderCells(combined_subset_pseudotime, root_state = root_state, num_paths = NULL, reverse = reverse) 
  
  ## re-do plotting based on refined root state
  cat ("re-do pseudotime trajectory plotting ...\n")
  pdf(paste0(out_path,"Updated_Pseudotime_along_trajectory_qval_",q_val_cutoff,".pdf"),width=8,height=8)  
  p <- plot_cell_trajectory(combined_subset_pseudotime, color_by = "Pseudotime",cell_size=0.3,cell_link_size=0.3,state_number_size=1)+theme(aspect.ratio=1)
  grid.draw(p)
  dev.off()
  
  ## save updated monocle object
  cat ("saving updated monocle object.\n")
  saveRDS(combined_subset_pseudotime,obj_output_path)
  
  sink()
} else {
  # Generate Monocle Object if It Does Not Exist Yet--------------------------------------------------
  ## Prepare data to Construct monocle object
  ### input merged/integrated object
  cat ("Read in Seurat object....\n")
  combined <- readRDS(path_srat)
  ### input annotation files
  idmetadata_df <- read.table(file = path_idmetadata, head=TRUE,sep="\t")
  barcode_to_celltype <- read.table(file = path_bc2celltype,head=TRUE,sep="\t")
  ### add annotation to the meta.data of the object
  #### get barcode in the combined object
  BC <- combined@meta.data %>% rownames
  #### get original barcode
  combined@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
  combined@meta.data <- merge(combined@meta.data %>% rownames_to_column, idmetadata_df[,c("Case", "Aliquot.snRNA","Sample_Type","Aliquot.snRNA.WU")],by.x="orig.ident",by.y="Aliquot.snRNA",all.x=TRUE)
  combined@meta.data <- merge(combined@meta.data,barcode_to_celltype,
                              by.x=c("orig.ident","original_barcode"), by.y=c("orig.ident","individual_barcode"))
  combined@meta.data <- column_to_rownames(combined@meta.data,"rowname")[BC,]
  ### subset the data to extract proximal tubules only
  cat ("Subset cells of interest....\n")
  combined_subset <- subset(combined,cells=combined@meta.data %>% rownames_to_column %>% filter(Cell_group.shorter=="Nephron_Epithelium") %>% .$rowname) 
  ### get gene count matrix
  exprs_matrix <- GetAssayData(combined_subset,slot="counts")
  ### get feature names
  fd <- data.frame(gene_short_name=rownames(exprs_matrix))
  rownames(fd) <- fd[,1]
  ### get phenotypic data
  pd <- combined_subset@meta.data
  
  ## Construct monocle object
  cat ("Construct monocle object (CellDataSet)....\n")
  combined_subset_pseudotime <- newCellDataSet(as(exprs_matrix,"sparseMatrix"),
                                               phenoData = new("AnnotatedDataFrame",data = pd),
                                               featureData = new("AnnotatedDataFrame",data = fd),
                                               lowerDetectionLimit = 0.5,
                                               expressionFamily=negbinomial.size())
  ## Estimate size factors and dispersions
  cat ("Estimate size factors...\n")
  combined_subset_pseudotime <- estimateSizeFactors(combined_subset_pseudotime)
  combined_subset_pseudotime <- estimateDispersions(combined_subset_pseudotime)
  ## Filter cells
  cat ("Filter Cells...\n")
  combined_subset_pseudotime <- detectGenes(combined_subset_pseudotime, min_expr = 0.1)
  print(head(fData(combined_subset_pseudotime)))
  expressed_genes <- row.names(subset(fData(combined_subset_pseudotime),num_cells_expressed >= 10))
  ## Clustering cells without marker genes
  cat ("Clustering..\n")
  disp_table <- dispersionTable(combined_subset_pseudotime)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  combined_subset_pseudotime <- setOrderingFilter(combined_subset_pseudotime, unsup_clustering_genes$gene_id)
  pdf(paste0(out_path,"ordering_genes_visualization.pdf"))
  p <- plot_ordering_genes(combined_subset_pseudotime)+theme(aspect.ratio=1)
  print(p)
  dev.off()
  
  combined_subset_pseudotime <- reduceDimension(combined_subset_pseudotime, max_components = 2, num_dim = 6,
                                                reduction_method = 'tSNE', verbose = T)
  combined_subset_pseudotime <- clusterCells(combined_subset_pseudotime, num_clusters = 6)
  
  
  
  pdf(paste0(out_path,"DimensionReductionPlot.pdf"))
  p1 <- plot_cell_clusters(combined_subset_pseudotime, 1, 2, color="Aliquot.snRNA.WU")+theme(aspect.ratio=1)
  #        p2 <- plot_cell_clusters(combined_subset_pseudotime, 1, 2, markers = c("Akap12","Sema5a","Il34","Cd44","Ccl2","Itgav"),cell_size=0.3)+theme(aspect.ratio=1)
  p3 <- plot_cell_clusters(combined_subset_pseudotime, 1, 2, color="seurat_clusters")+theme(aspect.ratio=1)
  print(p1)
  #        print(p2)
  print(p3)
  dev.off()
  
  ## Pseudotime trajectory construction
  ### select features for ordering
  cat ("Finding differentially expressed genes...\n")
  diff_test_res <- differentialGeneTest(combined_subset_pseudotime[expressed_genes,],
                                        fullModelFormulaStr = paste0("~",formula_str))
  write.table(diff_test_res,paste0(out_path,"diff_test_res.txt"),row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)
  
  ### use different qval as the thrshold for choosing ordering genes
  ordering_genes <- row.names (subset(diff_test_res, qval < q_val_cutoff))
  cat (paste0(length(ordering_genes)," features are used for ordering\n"))
  
  ## process according to ordering gene threshold (qval) as 0.01
  cat ("Ordering Cells...\n")
  combined_subset_pseudotime <- order_cells(monocle_obj=combined_subset_pseudotime,ordering_genes=ordering_genes,visualization_output=paste0(out_path,"ordering_genes_visualization2_qval",q_val_cutoff,".pdf"),reverse=FALSE)
  
  pdf(paste0(out_path,"plot_cell_Pseudotime_along_trajectory_qval_",q_val_cutoff,".pdf"),width=8,height=8)
  p <- plot_cell_trajectory(combined_subset_pseudotime, color_by = "Pseudotime",cell_size=0.3,cell_link_size=0.3,state_number_size=1)+theme(aspect.ratio=1)
  grid.draw(p)
  dev.off()
  
  saveRDS(combined_subset_pseudotime,obj_output_path)
  
  sink()
}

#pdf(paste0(out_path,"plot_cell_distribution_along_trajectory_qval_",q_val_cutoff,".pdf"),width=12,height=10)
#	plot_cell_trajectory(combined_subset_pseudotime, color_by = "Aliquot.snRNA.WU",cell_size=0.3,cell_link_size=0.3,state_number_size=1)+facet_wrap(~State, nrow = 3)+theme(aspect.ratio=1)
#	plot_cell_trajectory(combined_subset_pseudotime, color_by = "is_malignant",cell_size=0.3)+facet_wrap(~State, nrow = 1)+theme(aspect.ratio=1)
#	plot_cell_trajectory(combined_subset_pseudotime, color_by = "State",cell_size=0.3)+theme(aspect.ratio=1)
#dev.off()







