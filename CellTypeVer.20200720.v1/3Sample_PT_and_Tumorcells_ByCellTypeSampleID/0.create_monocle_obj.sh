# use cell type for differntial expression analysis; q cut-off set as 1e-10
id_run="3Sample_PT_and_Tumorcells_ByCellTypeSampleID"
path_out="/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Monocle/outputs/"${id_run}"/"
mkdir -p ${path_out}
Rscript ./pseudotime_analysis.R \
--FormulaStr=Cell_type.shorter+Aliquot.snRNA.WU \
--q_val=1e-10 \
--out_path=${path_out} \
--force.reprocess=TRUE &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
