# use cell type for differntial expression analysis; q cut-off set as 1e-10
source project_config.sh
mkdir -p ${path_out}
Rscript ./pseudotime_analysis.R \
--FormulaStr=Cell_type.detailed \
--q_val=1e-10 \
--out_path=${path_out} \
--force.reprocess=TRUE &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
