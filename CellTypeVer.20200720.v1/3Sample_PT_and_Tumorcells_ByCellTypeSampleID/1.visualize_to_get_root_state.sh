# use sampel id for differntial expression analysis; q cut-off set as 1e-10
id_run="3Sample_PT_and_Tumorcells_ByCellTypeSampleID"
path_out="/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Monocle/outputs/"${id_run}"/"
Rscript ./visualize_cell_distribution.R \
--out_path=${path_out} \
--monocle_obj=${path_out}combined_subset_pseudotime_qval_1e-10.rds &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
