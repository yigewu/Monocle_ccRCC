##-re-run step 0 with root state
source project_config.sh
Rscript ./pseudotime_analysis.R \
--FormulaStr=Cell_type.detailed \
--q_val=1e-10 \
--out_path=${path_out} \
--force.reprocess=FALSE \
--root_state=1 &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
