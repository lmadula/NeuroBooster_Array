#!/bin/bash

# Set for error handling and debugging

set -o pipefail   # Don't hide errors within pipes
set -x 		  # Enables bebugging

# Load the bcftools

spack load bcftools

# Define directories, PLINK path, QC thresholds, solved/unsolved cases, and annotation databases

Data_dir="$HOME/GP2_NBA_fixed"              # Directory containing input data
Output_dir="$HOME/fourth_qcoutput"          # Directory to store output files
Plink_path="$HOME/bin"                      # Path to the directory containing Plink
SNP_missingness="0.05"                      # SNP misssingness threshold
Individual_missingness="0.1"                # Individual missingness threshold 
Selected_vcf="$HOME/Unsolved_cases"         # Directory for unsolved cases VCF files
Combined_file="$Selected_vcf/combined.vcf"  # Merged VCF file for annotation
Annovar_db="$HOME/annovar/humandb"          # Path to ANNOVAR human database

# Base file name without extension

Base_name="flipped_GP2_merge_STELLENBOS"

# Input file paths

Bim_file="$Data_dir/${Base_name}.bim"
Fam_file="$Data_dir/${Base_name}.fam"
Master_key_file="$Data_dir/STELLENBOS_key_aug7_2024.csv"
PD_genes_file="$Data_dir/pd_gene_list.txt" 
Fasta_file="$Data_dir/Homo_sapiens_assembly38.fasta"

# Function to check if a file exists and has the correct extension

check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: File '$1' not found!" >&2
        exit 1
    fi

    # Check the file extension
    case "$1" in
        *.bim) ;;
        *.fam) ;;
        *.csv) ;;
        *.txt) ;;
        *.fasta) ;;
        *) 
            echo "Error: File '$1' has an invalid extension!" >&2
            exit 1
            ;;
    esac
}

# Check each file

check_file "$Bim_file"
check_file "$Fam_file"
check_file "$Master_key_file"
check_file "$PD_genes_file"
check_file "$Fasta_file"

# Proceed with processing if all files are valid

echo "All files are valid. Proceeding with processing..."


# Check if variants are present in the files

for var in '40310434' '40340400' '161350208' '161785820' '161350203' '20639990'; do
    grep "$var" "$Bim_file"
done

# Use masterkey to add phenotype and sex information for binary files 

awk -F',' '{print $1,$4}' "$Master_key_file" > "$Output_dir/pheno_info_mono.txt"

awk 'BEGIN {FS=OFS=" "}
    NR==FNR { phenotype[$1] = $2; next }
    { $6 = phenotype[$2]; print }' "$Output_dir/pheno_info_mono.txt" "$Fam_file" > "$Output_dir/updated.fam"

mv "$Output_dir/updated.fam" "$Output_dir/${Base_name}.fam"


awk -F',' '{print $1,$7}' "$Master_key_file" > "$Output_dir/sex_info_mono.txt"

awk 'BEGIN {FS=OFS=" "}
    NR==FNR { phenotype[$1] = $2; next }
    { $5 = phenotype[$2]; print }' "$Output_dir/sex_info_mono.txt" "$Fam_file" > "$Output_dir/updated.fam"

mv "$Output_dir/updated.fam" "$Output_dir/${Base_name}.fam"


# SNP Missingness 5%

$Plink_path/plink --bfile "$Data_dir/$Base_name" --geno "$SNP_missingness" --make-bed --out "$Output_dir/geno_$SNP_missingness"

 
# Check if variants are present in files after SNP fltering

for var in '40310434' '40340400' '161350208' '161785820' '161350203' '20639990'; do
    grep "$var" "$Output_dir/geno_$SNP_missingness.bim"
done

# Individual missingness 10%

$Plink_path/plink --bfile "$Output_dir/geno_$SNP_missingness" --mind "$Individual_missingness" --make-bed --out "$Output_dir/mind_$Individual_missingness"

# Check for sex discrepancies

$Plink_path/plink --bfile "$Data_dir/$Base_name" --check-sex 0.5 0.8 --out "$Output_dir/plink"
grep "PROBLEM" "$Output_dir/plink.sexcheck" | awk '{print$1,$2}' > "$Output_dir/sex_discrepancy.txt"
$Plink_path/plink --bfile "$Output_dir/mind_$Individual_missingness" --remove "$Output_dir/sex_discrepancy.txt" --make-bed --out "$Output_dir/sex_check_qc"

# Delete SNPs not in HWE
# Controls

$Plink_path/plink --bfile "$Output_dir/sex_check_qc" --hwe 1e-6 --make-bed --out "$Output_dir/hwe_controls"

# Cases

$Plink_path/plink --bfile "$Output_dir/hwe_controls" --hwe 1e-10 --hwe-all --make-bed --out "$Output_dir/hwe"

# Check if variants are present in files after HWE filtering
 
for var in '40310434' '40340400' '161350208' '161785820' '161350203' '20639990'; do
    grep "$var" "$Output_dir/hwe.bim"
done

# Runs of heterozygosity (LD)

$Plink_path/plink --bfile "$Output_dir/hwe" --indep-pairwise 50 5 0.2 --out "$Output_dir/indepSNP"
$Plink_path/plink --bfile "$Output_dir/hwe" --extract "$Output_dir/indepSNP.prune.in" --het --out "$Output_dir/R_check"

# Plot heterozygosity rate distribution

Rscript --vanilla -e "het <- read.table(file.path('$Output_dir', 'R_check.het'), header=TRUE);
pdf(file.path('$Output_dir', 'heterozygosity.pdf'));
het\$HET_RATE = (het\$'N.NM.' - het\$'O.HOM.') / het\$'N.NM.';
hist(het\$HET_RATE, xlab='Heterozygosity Rate', ylab='Frequency', main='Heterozygosity Rate');
dev.off()"

# Identify individuals deviating from the mean

Rscript --vanilla -e "het <- read.table(file.path('$Output_dir', 'R_check.het'), header=TRUE);
het\$HET_RATE = (het\$'N.NM.' - het\$'O.HOM.') / het\$'N.NM.';
het_fail = subset(het, (het\$HET_RATE < mean(het\$HET_RATE) - 3 * sd(het\$HET_RATE)) | (het\$HET_RATE > mean(het\$HET_RATE) + 3 * sd(het\$HET_RATE)));
het_fail\$HET_DST = (het_fail\$HET_RATE - mean(het\$HET_RATE)) / sd(het\$HET_RATE);
write.table(het_fail, file.path('$Output_dir', 'fail-het-qc.txt'), row.names=FALSE)"

# Adapt fail-het-qc.txt for PLINK

sed 's/"//g' "$Output_dir/fail-het-qc.txt" | awk '{print $1, $2}' > "$Output_dir/het_fail_ind.txt"

# Remove heterozygosity rate outliers

$Plink_path/plink --bfile "$Output_dir/hwe" --remove "$Output_dir/het_fail_ind.txt" --make-bed --out "$Output_dir/het"

# Python script to fix alleles

python3 - <<END
import os

# Input filename for the PLINK output (Heterozygosity)

dataIn = "het"

# PLINK command to convert the binary file to VCF format
# The --double-id option is used to treat the entire sample ID as a single identifier

os.system(f"$Plink_path/plink --bfile '$Output_dir/{dataIn}' --recode vcf-iid bgz --out '$Output_dir/fl' --output-chr chrMT --double-id")

# Index the generated VCF file

os.system(f"bcftools index '$Output_dir/fl.vcf.gz'")

# Fix reference alleles using bcftools
# The -Oz option compresses the output and -o specifies the output file name

os.system(f"bcftools +fixref '$Output_dir/fl.vcf.gz' -Oz -o '$Output_dir/fl_1.vcf.gz' -- -f '$Fasta_file' -m flip -d")

# Index the fixed VCF file

os.system(f"bcftools index '$Output_dir/fl_1.vcf.gz' -f")

# Convert the fixed VCF file back to PLINK binary format
# Again using --double-id to handle sample IDs correctly

os.system(f"$Plink_path/plink --vcf '$Output_dir/fl_1.vcf.gz' --make-bed --out '$Output_dir/fl_binary' --double-id")
END

# Check if variants are present in files after fixing reference

for var in '40310434' '40340400' '161350208' '161785820' '161350203' '20639990'; do
    grep "$var" "$Output_dir/fl_binary.bim"
done

# Extract PD genes

$Plink_path/plink --bfile "$Output_dir/fl_binary" --extract range "$PD_genes_file" --make-bed --out "$Output_dir/extract_genes"

# Check if variants are present in extracted gene files

for var in '40310434' '40340400' '161350208' '161785820' '161350203' '20639990'; do
    grep "$var" "$Output_dir/extract_genes.bim"
done

# Convert to .vcf

$Plink_path/plink --bfile "$Output_dir/extract_genes" --recode vcf --out "$Output_dir/extract_genes"


# Split samples in .vcf file

for sample in $(bcftools query -l "$Output_dir/extract_genes.vcf"); do
    bcftools view -c1 -s ${sample} "$Output_dir/extract_genes.vcf" -o "$Output_dir/${sample}.vcf"
done

# Select the unsolved cases    
# Define the pattern to match VCF files from 000001 to 000655

PATTERNS="STELLENBOS_000(00[1-9]|0[1-9][0-9]|[1-5][0-9]{2}|600|601|602|603|604|605|606|607|608|609|610|611|612|613|614|615|616|617|618|619|620|621|622|623|624|625|626|627|628|629|630|631|632|633|634|635|636|637|638|639|640|641|642|643|644|645|646|647|648|649|650|651|652|653|654|655)_s1\.vcf"

# Find and copy selected VCF files

find "$Output_dir" -type f -name "*.vcf" | grep -E "$PATTERNS" | while read -r file; do
   echo "Found and copying file: $file"
   cp "$file" "$Selected_vcf/"
done

# Combine the selected VCF files into one file

echo "Combining files into $Combined_file"
cat "$Selected_vcf"/*.vcf > "$Combined_file"

# Convert the combined VCF file to .avinput format

perl "$HOME/annovar/convert2annovar.pl" -format vcf4 "$Combined_file" > "$Selected_vcf/combined.avinput"

# Annotate using ANNOVAR

perl "$HOME/annovar/table_annovar.pl" "$Selected_vcf/combined.avinput" "$Annovar_db" -buildver hg38 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,clinvar,gnomad30_genome -operation g,r,f,f,f,f,f -nastring NA -csvout -polish

# Move the final annotation output to the selected VCF directory

mv "myanno.hg38_multianno.csv" "$Selected_vcf/myanno.hg38_multianno.csv"

## Analysis using R programming language

# Execute R script to prioritize variants


Rscript --vanilla -e "

# Load required libraries

library(dplyr)    # For data manipulation
library(readr)    # For reading CSV files
library(stringr)  # For string operations

# Load the dataset from ANNOVAR

myanno_hg38_multianno <- read_csv('${Selected_vcf}/myanno.hg38_multianno.csv', show_col_types = FALSE)

# Select essential variables

annovar_positive_controls <- myanno_hg38_multianno %>% 
                              select(Chr:Gene.refGene, AAChange.refGene, ExonicFunc.refGene,  
                                avsnp147, SIFT_score, Polyphen2_HDIV_score, 
                                  FATHMM_score, CADD_phred, CLNDN, -End, CLNSIG,AF,AF_afr,AF_eas,AF_nfe,AF_sas, MetaLR_score, MetaSVM_score)

# Renaming the second variable to be the same as the VCF variable (POS)

annovar_positive_controls <- annovar_positive_controls %>% rename(POS=Start)

# Filter rows containing 'Parkinson' or 'Parkinsonian'

Phenotype_PD <- annovar_positive_controls %>%
                 filter(str_detect(CLNDN, regex('Parkinson|Parkinsonian', ignore_case = TRUE))) %>%                          
                   mutate(Parkinson_disease = sapply(str_extract_all(CLNDN, 
                    regex('([^|]*Parkinson[^|]*|[^|]*Parkinsonian[^|]*)',
                    ignore_case = TRUE)), function(x) paste(x, collapse = ', ')))

# Define Filtering Criteria

# Define Filtering Criteria

metasvm_cutoff <- 0.8
metalr_cutoff <-  0.8

# Filter Variants based on MetaSVM OR MetaLR scores
 filtered_variants <- Phenotype_PD %>%
                      filter(MetaSVM_score > metasvm_cutoff | MetaLR_score > metalr_cutoff)


#cadd_cutoff <- 20
#polyphen_cutoff <- 0.5
#sift_cutoff <- 0.05 
#fathmm_cutoff <- 0

# Filter Variants based on in silico pathogenicity prediction tools

#filtered_variants <- Phenotype_PD %>%
 #filter(CADD_phred >= cadd_cutoff |
  # Polyphen2_HDIV_score > polyphen_cutoff |
   # SIFT_score < sift_cutoff |
    # FATHMM_score < fathmm_cutoff) %>% distinct()

	       
# Save the prioritized variants

write.csv(filtered_variants, file.path('${Selected_vcf}', 'prioritized_variants.csv'), row.names = FALSE)

# Define the path for the prioritized variants file

prioritized_variants_file <- file.path('${Selected_vcf}', 'prioritized_variants.csv')

# Read the prioritized variants file

prioritized_variants <- read.csv(prioritized_variants_file)

# Processing multiple VCF files

vcf_files <- list.files('${Selected_vcf}', pattern ='_s1*.vcf', full.names = TRUE)

# Loop through each VCF file and merge with the prioritized variants

all_merged_data <- data.frame()  # Initialize an empty data frame to store merged results

for (vcf_file in vcf_files) {
  # Read the current VCF file
  vcf_data <- read_delim(vcf_file, delim = '\t', escape_double = FALSE,
                         col_names = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GENOTYPE', 'HEROZYGOSITY_HOMOZYGOSITY'),
                         comment = '#', trim_ws = TRUE)

  # Add a column for the VCF source dynamically using the file name (basename)

  vcf_data <- vcf_data %>%
               mutate(VCF_SOURCE = basename(vcf_file))

  # Ensure position is numeric in both datasets

  prioritized_variants\$POS <- as.numeric(as.character(prioritized_variants\$POS))
  vcf_data\$POS <- as.numeric(as.character(vcf_data\$POS))

  # Remove rows with NA values in position from both datasets

  prioritized_variants <- prioritized_variants %>% filter(!is.na(POS))
  vcf_data <- vcf_data %>% filter(!is.na(POS))

  # Filter the prioritized variants based on matching positions with the VCF

  aligned_data <- prioritized_variants %>%
                    filter(POS %in% vcf_data\$POS)

  # Merge the filtered CSV data with the current VCF data based on position

  merged_data <- aligned_data %>%
                  left_join(vcf_data, by = 'POS')

  # Append the merged data to the all_merged_data data frame

  all_merged_data <- bind_rows(all_merged_data, merged_data)
}

# Save final merged data

write.csv(all_merged_data, file.path('${Selected_vcf}', 'merged_data.csv'), row.names = FALSE)"









