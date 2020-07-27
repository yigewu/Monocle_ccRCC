# use sampel id for differntial expression analysis; q cut-off set as 1e-10
path_out=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Monocle/outputs/C3L-00088_test/
Rscript ./pseudotime_analysis_for_C3L-00088_test.R \
--FormulaStr=Aliquot.snRNA.WU \
--q_val=1e-10 \
--out_path=${path_out} \
--force.reprocess=TRUE &> ${path_out}Log.$(date +%Y%m%d%H%M%S).txt&
