# APAP_toxicity_liver
# Genome-scale model analysis of acetaminophen toxicity in the laboratory rat 
Code availability for the manuscript titled "Metabolic network-based predictions of toxicant-induced metabolite changes in the laboratory rat" by Pannala *et al*.

Following script was used to generate Figures in the main text and Supplementary Fig. S1. 

## Figure 2
Inputs processed gene expression fold changes stored in the folder Weights_apap
Outputs volcano plot for gene expression changes

## Figure 4
Inputs processed metabolic flux analysis data
Outputs bar charts for absolute flux changes

## Figure 5
Inputs processed metabolite fold changes reported in Suppplementary TableS3
Outputs volcano plots for metabolite changes

## Figure 7
Inputs processed table of TIMBR predcitions for the metabolite changes given in the folders MFA_10h and NoMFA_10h and observed metabolite changes Mets_increased/decreased_10h_data
Outputs heatmap showing the observed metabolite changes in the data vs. model predictions at 10 h post APAP treatment
Helper functions: Timbr_10h_MFA.m and Timbr_10h_NoMFA.m are used to generate the Tables in the MFA_10h folders

## Figure S1
Inputs processed table of TIMBR predcitions for the metabolite changes given in the folders MFA_5h and NoMFA_5h and observed metabolite changes Mets_increased/decreased_10h_data
Outputs heatmap showing the observed metabolite changes in the data vs. model predictions at 5 h post APAP treatment 
Helper functions: Timbr_5h_MFA.m and Timbr_5h_NoMFA.m are used to generate the Tables in the MFA_5h folders

## Predicting biomarkers based on gene expression changes with __TIMBR__

### Generate_timbr_weights.R

  Inputs gene expression changes from the folder Weights_apap, Sleuth_table_liver_.tsv files, Supplementary_iRno_v2.xlsx and gene annotation tables from the file, Biomart_rno_ids.txt.gz
  
  Outputs TIMBR reaction weights
  
  Helper functions: ncomm_helper.R

### Timbr_h_.m

  Inputs TIMBR reaction weights, rat metabolic network in the cobra format iRno_v2.mat
  
  Outputs TIMBR predictions for metabolite changes
  
  Helper functions: timbr.m and ncomm_blais_model2irrev.m

## iRno_v2.mat
  Rat genome scale metabolic model version 2 in cobra format

## iRno_v2_SBML.xml
  Rat genome scale metabolic model version 2 in SBML

## Helper function files:

### ncomm_helper.R
  R source file that defines several helper functions used in the R scripts above.

### timbr.m
  MATLAB implementation of the TIMBR (transcriptionally-inferred metabolic biomarker response) algorithm

### ncomm_blais_model2irrev.m
  Modified version of convertToIrreversible from the COBRA toolbox (www.github.com/opencobra/cobratoolbox) 
  that facilitates the mapping of TIMBR reaction weights to irreversible reactions
