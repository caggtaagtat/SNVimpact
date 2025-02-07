# 5'ss variations

This repository contains code to reproduce the analysis of our paper with the titel: 
### Fully haplotyped genome assemblies of healthy individuals reveal variability in 5’ss strength and support by splicing regulatory proteins 

Supplementary Data and data used for the analysis can be accessed via zonodo und DOI: X

First, the large VCF file from the publication was preprocessed into 2 csv tables stating haplotyped sequence variations per individual genome (1_Preprocessing VCF.sh, 2_Create_csv.R)

Next, every identifed sequence variation was introduced into the reference sequence surrounding of 5'ss (3_Introduce_variations.sh), applying 3 custom R scripts (31_Introduce_variations.R, 32_SDvariation.R and 33_GetReferenceSequence.R).

The information of the indiviudal csv-files was then summarized and analyzed with an additional custom R-script (4_Evaluation_5ss_variation.R). All figures of the paper can be generated with the final dataset (Final_data.R) and the custom R-script (5_Figure_generation.R).
