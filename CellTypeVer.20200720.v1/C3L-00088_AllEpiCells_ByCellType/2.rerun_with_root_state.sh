##-re-run step 0 with root state
path_out=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Monocle/outputs/C3L-00088/
Rscript ./pseudotime_analysis.R \
--FormulaStr=Cell_type.shorter \
--q_val=1e-10 \
--out_path=${path_out} \
--force.reprocess=FALSE \
--root_state=5 &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
