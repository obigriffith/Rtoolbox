########################################################################
# Prediction of response to 90 targeted and cytotoxic agents based on	 #
# gene expression data, methylation data and/or copy number data       #
#                                                                      #
# Copyright 2012 Anneleen Daemen, Obi Griffith, Laura Heiser           #
########################################################################

We provide an R script through which the list of 90 drug compounds is ranked for a patient according to predicted drug response based on gene expression, copy number and/or methylation data, upon availability for that patient. Best results are expected if data are provided from platforms used for the cell line data, based on which the classifiers were built (that is, Affymetrix's GeneChip Human Genome U133A and Genome-Wide Human SNP Array 6.0 for expression and copy number data, respectively, and Illumina's Infinium HumanMethylation BeadChip for methylation). In that case, use R script DrugResponseToolboxPlatformSpecific.R. The input files per patient should be one or more of the following: U133A CEL-file, SNP6 CEL-file, and a tab-delimited file with the proportion of methylated DNA at each measured CpG locus. The U133A CEL-file is preprocessed with the U133A CEL-files of the 48 core cell lines, whilst the SNP6 CEL-file is segmented using the same breakpoints as obtained in the cell line panel. Predicted response for each patient is provided for those compounds with AUC>0.6 for the best cell line-derived predictor based on any or all of the input data.

However, we also provide a platform agnostic toolbox in R that is independent of the platforms used to obtain expression, copy number and/or methylation data. In this case, use R script DrugResponseToolboxAgnostic.R. The input files per patient should be one or more of the following: a tab-delimited file with gene symbol and expression level, a tab-delimited file with gene symbol and copy number level, and a tab-delimited file with the proportion of methylated DNA at each measured CpG locus. Both the expression and copy number data are first quantile normalized to the cell line distributions. Predicted response for each patient is provided for those compounds with AUC>0.6 for the best cell line-derived predictor based on any or all of the input data and with provided patient data on a sufficient number of predictor/model variables. The latter requires a weighted percent of model variables (WPMV) of at least 80% calculated as the sum of input (user-provided) variable ranks (where most important variables receive the largest rank value) divided by the sum of all model variable ranks.

NOTE: This github version of the source code is provided to make a version available without any special login requirements. However, we strongly recommend you do create an account and try the software from the Synapse repository:
https://www.synapse.org/#!Synapse:syn2179898


1) Installations required in R before running the R script (preferably with R version 2.12.0 or 2.12.2):

source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("affy")
biocLite("DNAcopy")
biocLite("CNTools")
biocLite("preprocessCore")

source("http://aroma-project.org/hbLite.R")
hbInstall("aroma.affymetrix")
hbInstall("aroma.cn")

install.packages(randomForest)
install.packages(gdata)


2) Verify whether pdflatex is available on your system

If not, MacTex can be installed from http://www.tug.org/mactex/2010/. The smaller package BasicTeX.pkg.zip is sufficient. For Centos linux, pdflatex is packaged in tetex-latex. In Ubuntu linux it is part of tetex-extra.


3) Folders that are provided in this package

- Folder CoreCELfiles: CEL files of the core cell lines, used for the preprocessing of each individual tumor CEL file
- Folder ExtraFiles: contains additional files required for preprocessing and model testing - list of drug compounds, list of best performing models per drug depending on the input files, methylation cell line data and annotation, U133A cell line data, and a SNP6 subfolder with files required for the preprocessing
- Folder ModelsPlatformSpecific and ModelsAgnostic: Best predictors for the 90 drug compounds, with predictor based on expression, methylation and/or copy number data, depending on the patient's input files
- Folder TumorData: The data files for the patient(s) should be saved to this folder
- Folder Results: results will be saved to this folder when running the R script

NOTE: Because of github file size limitations some files will have to be downloaded separately.
Make sure that you are in the main 'Rtoolbox' folder and then run the following commands:
wget https://dl.dropboxusercontent.com/u/16769159/BCCL/CoreCELfiles.tar.gz
wget https://dl.dropboxusercontent.com/u/16769159/BCCL/ExtraFiles.tar.gz
tar -zxvf CoreCELfiles.tar.gz
tar -zxvf ExtraFiles.tar.gz

4) Changes to be made in the R script DrugResponseToolbox.R

In the beginning of the R script, path needs to be changed to the user's specific path.
* path = path of the unzipped package, containing the R script and the folders listed under 3\
The remaining of the script is relative to this path.


5) Running the R script

The R script should be called on command line with as arguments the patient ID and the names of the tumor CEL or txt file(s), with "NA" in case of unavailability of one or two particular data types with order U133A - methylation - SNP6, for example:
> ./DrugResponseToolboxPlatformSpecific.R "Cellline" "U133A_test_tumor.CEL" "NA" "SNP6_test_cellline.CEL"
> ./DrugResponseToolboxAgnostic.R "Patient" "Expression_test_tumor.txt" "Meth_test_tumor.txt" "CNV_test_tumor.txt"


6) Result-files generated per tumor sample

- txt-file per patient with the table of scores, probability of response, sensitivity, and data combination/model the prediction was based on for the 90 drug compounds
- pdf-file per patient, containing the names of the 1-3 input files, the same table as above and recommended treatment
