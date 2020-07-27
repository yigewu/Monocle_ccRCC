# use sampel id for differntial expression analysis; q cut-off set as 1e-10
path_out=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Monocle/outputs/C3L-00088/
Rscript ./visualize_cell_distribution.R \
--out_path=${path_out} \
--monocle_obj=${path_out}combined_subset_pseudotime_qval_1e-10.rds &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
