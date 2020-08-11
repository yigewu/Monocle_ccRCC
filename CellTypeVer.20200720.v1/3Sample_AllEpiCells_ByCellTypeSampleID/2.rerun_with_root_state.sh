##-re-run step 0 with root state
id_run="3Sample_AllEpiCells_ByCellTypeSampleID"
path_out="/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Monocle/outputs/"${id_run}"/"
Rscript ./pseudotime_analysis.R \
--FormulaStr=Cell_type.detailed+Aliquot.snRNA.WU \
--q_val=1e-10 \
--out_path=${path_out} \
--force.reprocess=FALSE \
--root_state=1 &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
