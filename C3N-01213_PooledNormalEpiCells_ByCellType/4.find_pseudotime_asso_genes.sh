## !! DO NOT RUN THIS IF THE MONOCLE OBJECT IS NOT UPDATED; TAKES A HUGE AMOUNT OF MEMORY AND TIME
source project_config.sh
Rscript ./finding_pseudotime_associated_genes.R \
--branch_point=1,2,3 \
--out_path=${path_out}/pseudotime_associated_genes/ \
--monocle_obj=${path_out}combined_subset_pseudotime_qval_1e-10.rds &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
