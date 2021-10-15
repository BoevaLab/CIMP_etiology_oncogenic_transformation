# Deciphering the etiology and role in oncogenic transformation of the CpG island methylator phenotype (CIMP): a pan-cancer analysis
This repository contains the Jupyter Notebooks/scripts used for the analysis in "Deciphering the etiology and role in oncogenic transformation of the CpG island methylator phenotype (CIMP): a pan-cancer analysis".

### How to
The code presented is shared for reproducibility purposes and is in no way tailored for wide-spread use, nor is it organized as a package. Please contact josephine.yates@inf.ethz.ch for more information.

**IMPORTANT DISCLAIMER**: a number of files are used in the analysis - placeholders have been integrated in their place in the provided code. Note that to run the code one will have to download/generate all necessary files and change the placeholders. A thorough description of where to download associated files is provided beneath, however there might be some slight changes to the structure of the files etc. which might cause the processing not to run. As this is probably the result of a change in the structure of the file, it can be probably resolved by verifying confirmity of the structure with the code. If there are any issues, please contact josephine.yates@inf.ethz.ch.

#### Preprocessing
The notebook `Deciphering_CIMP_Preprocessing.ipynb` is used for the thorough preprocessing of the data (non-CGI, NA, non-variable, purity correction and age-related). More details can be found in the Methods section of the aforementioned paper. 

What data to download to replace the placeholders: 
- survival_dir = files can be downloaded from the UCSC Xena data portal, under the sections "phenotype - survival data" for all analyzed cancer types. For example, the file for ACC is [this file](https://xenabrowser.net/datapages/?dataset=TCGA-ACC.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). 
- clin_dir = files can be downloaded from the UCSC Xena data portal, under the sections "phenotype - phenotype" for all analyzed cancer types. For example, the file for ACC is [this file](https://xenabrowser.net/datapages/?dataset=TCGA-ACC.GDC_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).
- info_dir = path to the folder to save the results of the preprocessing of clinical and methylation data.
- raw_dir = files can be downloaded from the UCSC Xena data portal, under the sections "DNA methylation - Illumina Human Methylation 450" for all analyzed cancer types. For example, the file for ACC is [this file](https://xenabrowser.net/datapages/?dataset=TCGA-ACC.methylation450.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). 
- estimate_purity_dir = ESTIMATE purity estimations can be downloaded from [here](https://bioinformatics.mdanderson.org/estimate/disease.html). We want to save the ESTIMATE scores (preprocessing might vary depending on the format of the downloaded file).
- CPE purity file = the file corresponds to the "Supplementary Data 1" file of [this article](https://www.nature.com/articles/ncomms9971). The file needs to be saved in csv format prior to being uploaded in the notebook.
- Illumina 450k annotation file = the file corresponds to "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz" and can be downloaded from [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534). 
- std_dir, purity_dir, nnls_dir, age_r_dir = folders to save your analyses for standard deviation filtering, purity correction, NNLS estimates, age related cpgs
- admp_dir = the files in this directory correspond to the age-dependent cpgs identified by Slieker et al. in [their paper](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-018-0191-3#Sec15), that can be found in Table S2. For ease, this file was transfored in the following manner: each tissue-specific loss (resp. gain) CpG was saved in a .txt file where every CpG probe was separated by `\n` named `TISSUE_loss.txt` (resp. `_gain`).

#### Analysis
The notebook `Deciphering_CIMP_Analysis.ipynb` is used for the subsequent analysis as described in the Methods section.

What data to download to replace the placeholders: 
- info_dir = path to the folder where you saved your preprocessed clinical data.
- clin_raw_data = path to the folder where you saved the raw clinical data (cf. clin_dir in the previous paragraph).
- meth_dir = path to the folder where you saved your preprocessed methylation data.
- mut_dir = files can be downloaded from the UCSC Xena data portal, under the sections "somatic mutation (SNP and INDEL) - MC3 gene-level non-silent mutation" for all analyzed cancer types. For example, the file for ACC is [this file](https://xenabrowser.net/datapages/?dataset=mc3_gene_level%2FACC_mc3_gene_level.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). *IMPORTANT: these files are in the TCGA dataset and not the GDC dataset as the previous files were.*
- msi profile file = the file corresponds to the Supplementary Data 1 in [this paper](https://www.nature.com/articles/ncomms15180). This should be saved as a csv file.
- geo_normal_dir = these are the transformed files from GEO accession series [GSE77871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77871) and [GSE32149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32149). The preprocessed files are given in the repository in the zip file in the `geo_normal_tissue` folder.
- cluster_memb_dir, sign_cpg_dir, bac_pred_dir = path to the folder to save the cluster memberships, significantly differentially methylated probes, balanced accuracy predictions in
- gencode_annot_gene and merged_table = preprocessed gencode annotation mappings that give the associated length and the mapping between the old and new versions of gencode annotations. These are given in the repository in the zip file in the folder `gencode_annotations`. 
- immune_scores, CIBERSORT_comp = these files correspond to files generated in [this paper](https://www.sciencedirect.com/science/article/pii/S1074761318301213?via%3Dihub). The files were downloaded from the GDC [here](https://gdc.cancer.gov/about-data/publications/panimmune). Immune_scores corresponds to the `Scores_160_Signatures.tsv.gz` file. CIBERSORT_comp corresponds to the `TCGA.Kallisto.fullIDs.cibersort.relative.tsv` file. 
-  thorsson_pat_info = this file corresponds to the file generated in [this paper](https://www.sciencedirect.com/science/article/pii/S1074761318301213?via%3Dihub). It corresponds to the Table S1 file.
-  diff_gex_dir = this is the folder where the results of DESeq2 analysis are saved (see paragraph below).
-  diff_gex_offgene_dir = folder to save the official gene codes for genes that correspond to "potential downstream events" (cf. Methods of the paper). *IMPORTANT: these are the files that are inputted into the [DAVID](https://david.ncifcrf.gov/summary.jsp) tool for the gene set enrichment analysis.*
-  tss_dir = folder to save the TSS distances generated
-  gencode_annot_dir = folder where the gencode annotations are downloaded. These can be taken from [here](https://www.gencodegenes.org/human/). They correspond to the .gtf file.

#### Differential Gene Expression analysis
With the results of the clustering membership, we computed a Differential Gene Expression analysis. This analysis was conducted in R so as to be able to use the tool DESeq2. The script used to generate the differential gene expression analysis is `DESeq_GEX.R`. 

What data to download to replace the placeholders: 
- directory = files can be downloaded from the UCSC Xena data portal, under the sections "gene expression RNAseq - HTSeq - Counts" for all analyzed cancer types. For example, the file for ACC is [this file](https://xenabrowser.net/datapages/?dataset=TCGA-ACC.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).
- file_clust = the folder where the cluster memberships generated with `Deciphering_CIMP_Analysis.ipynb` are saved.

