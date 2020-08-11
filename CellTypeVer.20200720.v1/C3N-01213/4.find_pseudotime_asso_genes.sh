## !! DO NOT RUN THIS IF THE MONOCLE OBJECT IS NOT UPDATED; TAKES A HUGE AMOUNT OF MEMORY AND TIME
id_run="C3N-01213"
path_out="/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Monocle/outputs/"${id_run}"/"
Rscript ./finding_pseudotime_associated_genes.R \
--branch_point=1,2 \
--out_path=${path_out}/pseudotime_associated_genes/ \
--monocle_obj=${path_out}combined_subset_pseudotime_qval_1e-10.rds &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
