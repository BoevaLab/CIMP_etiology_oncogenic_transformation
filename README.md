# CIMP_etiology_oncogenic_transformation
This repository contains the Jupyter Notebooks used for the analysis in "Deciphering the etiology and role in oncogenic transformation of the CpG island methylator phenotype (CIMP): a pan-cancer analysis".

### How to
The code presented is shared for reproducibility purposes and is in no way tailored for wide-spread use, nor is it organized as a package. Please contact josephine.yates@inf.ethz.ch for more information.

The notebook `Deciphering_CIMP_Preprocessing.ipynb` is used for the thorough preprocessing of the data (non-CGI, NA, non-variable, purity correction and age-related). More details can be found in the Methods section of the aforementioned paper

The notebook `Deciphering_CIMP_Analysis.ipynb` is used for the subsequent analysis as described in the Methods section.

**IMPORTANT DISCLAIMER**: a number of files are used in the analysis - placeholders have been integrated in their place in the provided code. All the files used for analysis are either available directly online as described in the Methods section or can be derived from the online sources (*eg* purity estimates from ESTIMATE can  be derived from publically available TCGA data). Note that to run the code one will have to generate all necessary files and change the placeholders.
