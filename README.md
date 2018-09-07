# UnlabelledProteomicsImputation
Imputation of missing values for unlabelled MS proteomics data 


The repository currently contains Python 3 script "gen_pr_files_3.py" for initial preparation 
of proteomics and transcriptomics files and R script "imputation_models_M1_M2_M3_M4.R" for
imputing missing protein values using one of the regression models: M1, M2, M3 or M4.

Sample proteomics and transcriptomics files for NCI-60 data cell line data set are also included.


Python script "gen_pr_files_3.py": 

The purpose of the script is from given two files with proteomics and transriptomics data produce 
"well formatted" files that can be user by R imputation script - meaning that both files processed
by imputation script contain identical number of rows (with the i-th row in each file corresponding 
to the same gene) and identical number of columns (with the j-th column in each file corresponding 
to the same tissue type/cell line). 
To obtain useful output the input files both should contain: 1) column (identical header names in both
files) containing IDs of genes, 2) at least one column each (identical header names in both files)
with quantitative proteomics and transcriptomics data correspondingly, 3) at least one row each with
the same gene ID. The output files will contain the maximal number of columns and rows that will
be possible to match acccording to header names and gene ids with columns and rows sorted in the same order
in both files. In case of ambiguous input data (e.g. two columns with the same header names) the output
result will be chosen randomly.

A brief description of "gen_pr_files_3.py" is given in "gen_pr_files_3_ description.py".
The sample input files for this script are "nci60_proteomics_raw.tsv" and "nci60_transcriptomics_raw.tsv". 
The script on these input files can be run with command:

python gen_pr_files_3.py -i nci60_proteomics_raw.tsv -I nci60_transcriptomics_raw.tsv -o nci60_prot.tsv -O nci60_rna.tsv -e "Gene_id" -f "TAB" 1> logfile.txt 2> errorfile.txt        


R script "imputation_models_M1_M2_M3_M4.R":

The script require "well formatted" input files of transcriptomics and proteomics data (sample data files
provided are "nci60_prot.tsv" and "nci_rna.tsv". Non-numerical NA and 0 values all will be converted to 
0-es, and treated as "real 0-es" in transcriptomics file and as missing values to be imputed in proteomics 
file. The default output files are "out_correlations.tsv" showing performance for each of the models and 
"out_imputed.tsv" with non-numerical, NA and 0 values in proteome file replaced with the imputed ones. In
addition the predicted negative values will be replaced by 0-es. 
The default imputation model is M4 with 2 iterations. Input and output file namaes and imputation model can be 
changed by simple edits in initial header part of the script.

 

 