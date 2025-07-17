🧠 Computational Workflow for Parkinson’s Disease Screening
Developed during MSc Research
This repository contains the complete computational workflow I developed during my Master’s research project focused on screening for neurodegenerative disorders, with a specific emphasis on Parkinson’s disease.

The analysis is based on genotyping array data, and the workflow is designed to streamline the process of data quality control, annotation, and statistical analysis using widely adopted bioinformatics tools.

📁 Folder Structure & Scripts
Bash_script/
This folder contains the main Bash script that automates the core steps of the analysis workflow, including:
Quality control and filtering using PLINK
Data formatting and preprocessing
Variant annotation using ANNOVAR
Exporting the cleaned and annotated data into a format suitable for downstream statistical analysis in R

R/
This folder contains the R script titled MSc_analysis.R, which performs:
Data import 
Visualization of key results (e.g., Number of variants per gene/chromosome, amino acid change, allele frequency, etc.)
Summary tables and figures

📊 Data Description
The analysis is based on SNP array genotype data derived from case-control samples relevant to Parkinson’s disease. The dataset includes individual-level genotype calls, which were processed to identify variants potentially associated with disease phenotypes.

🎯 Objectives
The primary goals of this workflow are to:
Conduct a reliable and reproducible quality control of genotyping data
Functionally annotate variants 
Identify variants associated with Parkinson’s disease (By filtering based Parkinson/Parkinsonian and prioritised the variants using _in-silico_ prediction tools)
Provide a structured and automated workflow that can be reused or adapted for related studies

🔁 Reproducibility
All scripts are documented and organized to ensure reproducibility. To run the complete workflow:
Start with the Bash script located in Bash_script/
Follow the automated steps for quality control and annotation
Export the processed data and run the MSc_analysis.R script in the R/ folder for further analysis

📜 License
This project is released under the MIT License.




