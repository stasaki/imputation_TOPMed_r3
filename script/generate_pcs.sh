#!/bin/bash

cd /data/shinya/shinya/Project/imputation_TOPMed_r3
conda activate gwas_imputation
# Define file paths
input_prefix="./combined_plink/TOPMed_r3_all_hg38"
output_prefix="./combined_plink/TOPMed_r3_all_hg38"
maf_threshold=0.01
hwe_threshold=1e-6
geno_threshold=0.05
ld_window_size=50
ld_step_size=5
ld_r2_threshold=0.2

# Step 1: Quality Control
plink --bfile $input_prefix --maf $maf_threshold --hwe $hwe_threshold --geno $geno_threshold --make-bed --out ${output_prefix}_qc

# Step 2: Linkage Disequilibrium (LD) Pruning
plink --bfile ${output_prefix}_qc --indep-pairwise $ld_window_size $ld_step_size $ld_r2_threshold --out ${output_prefix}_prune

# Step 3: Generate Principal Components (PCs)
plink --bfile ${output_prefix}_qc --extract ${output_prefix}_prune.prune.in --pca --out ${output_prefix}_pca

echo "PCA analysis completed. PCs saved in ${output_prefix}_pca.eigenvec and ${output_prefix}_pca.eigenval."
rm ${output_prefix}_qc.*
rm ${output_prefix}_prune.*
rm ${output_prefix}_pca.log
rm ${output_prefix}_pca.nosex 
rm ${output_prefix}_pca.eigenval 
# run filter_samples.R
rm ${output_prefix}_pca.eigenvec