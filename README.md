# TCGA_BRCA_ERBB2_Expression_Analysis
Gene expression analysis of ERBB2-amplified breast cancer using TCGA data
TCGA BRCA ERBB2 Gene Expression Analysis

This repository contains the code and outputs used for Assignment 2: Gene Expression Analysis and Interpretation, focusing on ERBB2 (HER2) amplification in breast cancer using TCGA data.

Data used

Publicly available TCGA breast cancer (BRCA) data were downloaded from cBioPortal:

RNA-seq gene expression data (data_mrna_seq_v2_rsem.txt)

Copy number alteration (CNA) data (data_cna.txt)

Clinical patient data (data_clinical_patient.txt)

Code overview

The analysis was developed interactively in R and saved as a single script.

analysis.R (or console_history.R)

This script performs the complete analysis workflow:

Data loading and cleaning

Reads RNA-seq, CNA, and clinical data files

Removes metadata rows from clinical data

Ensures consistent patient identifiers across datasets

Patient matching and metadata creation

Matches patients across RNA-seq, CNA, and clinical datasets

Extracts ERBB2 copy number values

Defines ERBB2 amplification status (CNA > 0 = amplified)

Differential gene expression analysis

Prepares count data for DESeq2

Performs differential expression analysis between ERBB2-amplified and non-amplified tumours

Outputs full and filtered lists of differentially expressed genes

Exploratory analysis and visualisation

Applies variance stabilising transformation (VST)

Generates PCA plot

Generates heatmap of top differentially expressed genes

Survival analysis

Extracts overall survival time and status from clinical data

Performs LASSO-regularised Cox proportional hazards regression using glmnet

Computes risk scores and stratifies patients into high- and low-risk groups

Generates Kaplan–Meier survival plot

Output files
Output_figures/

Contains figures generated during the analysis:

PCA plot of variance-stabilised expression data

Heatmap of top differentially expressed genes

Kaplan–Meier survival curves for risk groups

Output_tables/

Contains tables generated during the analysis:

Matched metadata and expression matrices

Differential expression results (full and filtered)

Top differentially expressed genes

Genes selected by Cox LASSO regression

