#conda activate gwas_imputation
export PLINK_HOME=~/miniconda3/envs/gwas_imputation/bin/
 
cd /data/shinya/shinya/Project/imputation_TOPMed_r3/AMPAD_WGS/tmp
BEDNAME="NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08.recalibrated_variants_hg38_dn.QCed1"

cd /data/shinya/shinya/Project/imputation_TOPMed_r3/Diverse_Joint/tmp
BEDNAME="DIVERSE.joint.callset"

# LD pruning parameters
WINDOW_SIZE=200  # Size of the window in SNPs
STEP_SIZE=5     # Number of SNPs to shift the window at each step
R2_THRESHOLD=0.9  # R^2 threshold for pruning
MAF=0.01
VCF_HOME="/data/shinya/shinya/Project/imputation_TOPMed_r3/vcftools/src/perl/"
# Generate chromosome rename mappings only if needed
OUTPUT_FILE="chr_rename.txt"
if [ ! -f "$OUTPUT_FILE" ]; then
    > "$OUTPUT_FILE"
    for i in {1..22}; do
        echo "$i chr$i" >> "$OUTPUT_FILE"
    done
    echo "X chrX" >> "$OUTPUT_FILE"
else
    echo "Chromosome rename mappings already generated. Skipping."
fi


for i in `seq 1 22`; do
       if [ ! -f "../${BEDNAME}_filtered-updated-chr${i}.vcf.gz" ]; then
            # Define input and output files
            INPUT_DATA="${BEDNAME}_filtered-updated-chr${i}"  # Replace with your input file prefix (without extensions)
            OUTPUT_FILTERED="${BEDNAME}_filtered-updated-chr${i}-pruned"
            
            # Run PLINK to prune variants based on LD
            plink --bfile ${INPUT_DATA} --maf ${MAF} --indep-pairwise ${WINDOW_SIZE} ${STEP_SIZE} ${R2_THRESHOLD} --out pruned_data
            
            # Step 3: Create a new dataset with filtered SNPs
           #plink --bfile ${INPUT_DATA} --extract pruned_data.prune.in --make-bed --out ${OUTPUT_FILTERED}
            
            plink --bfile ${INPUT_DATA} --extract pruned_data.prune.in  --real-ref-alleles --recode vcf --out ${OUTPUT_FILTERED}

            rm pruned_data.prune.in
            rm pruned_data.prune.out
            rm pruned_data.hh
            rm pruned_data.log
            
            $PLINK_HOME/bcftools annotate --rename-chrs chr_rename.txt ${BEDNAME}_filtered-updated-chr${i}-pruned.vcf -Ov -o ${BEDNAME}_filtered-updated-chr${i}-pruned_renamed.vcf
            $VCF_HOME/vcf-sort ${BEDNAME}_filtered-updated-chr${i}-pruned_renamed.vcf | bgzip -c >  ../${BEDNAME}_filtered-updated-chr${i}.vcf.gz
            rm ${BEDNAME}_filtered-updated-chr${i}-pruned_renamed.vcf  
            rm ${BEDNAME}_filtered-updated-chr${i}-pruned.*  
        else
            echo "VCF file for chromosome ${i} already processed. Skipping."
        fi
done    






