#!/bin/bash
set -e
# Check if the correct number of arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bucket_location> <STUDY>"
    exit 1
fi

# Assign arguments to variables
BUCKET_LOCATION=$1
STUDY=$2

# Extract BEDNAME from BUCKET_LOCATION
BEDNAME=$(basename $BUCKET_LOCATION)

#conda activate plink

# Set environment variables
export PLINK_HOME=~/miniconda3/envs/plink/bin/
export IC_HOME=/home/shinya/Resource/GENETICS_Resource/softwares/imputation_check/pre/v4.3.0/
export REF=/home/shinya/Resource/GENETICS_Resource/softwares/imputation_check/ref/20240317/PASS.Variantsbravo-dbsnp-all.tab.gz
export checkVCF_HOME=~/Resource/GENETICS_Resource/softwares/checkVCF
export genome_LOC=~/Resource/NGS_Resource/genome/GRCh38.primary_assembly.genome.fa
export API_key=""

# Make a temporary directory and move into it
mkdir -p ${STUDY}/tmp
cd ${STUDY}/tmp

# Download data if not already present
if [ ! -f "${BEDNAME}.bim" ]; then
    gsutil cp ${BUCKET_LOCATION}.bim ./
    gsutil cp ${BUCKET_LOCATION}.fam ./
    gsutil cp ${BUCKET_LOCATION}.bed ./
else
    echo "${BEDNAME} files already downloaded. Skipping download."
fi



# Check for the presence of filtered output files to avoid re-running processes
if [ ! -f "${BEDNAME}_filtered.bim" ]; then
    # Pre-imputation quality control and file preparation
    #$PLINK_HOME/plink --bfile $BEDNAME --remove ../../do_not_distribute.txt --freq --make-bed --out ${BEDNAME}_filtered
    $PLINK_HOME/plink --bfile $BEDNAME --freq --make-bed --out ${BEDNAME}_filtered
    nohup perl $IC_HOME/HRC-1000G-check-bim.pl -b ${BEDNAME}_filtered.bim -f ${BEDNAME}_filtered.frq -r $REF -h -v
    nohup bash Run-plink.sh
else
    echo "Filtered files for ${BEDNAME} already exist. Skipping PLINK processing."
fi

# Proceed with checks only if necessary files are missing
if [ ! -f "${BEDNAME}_filtered-updated-chr22.check.log" ]; then
    # Check VCF files for errors and processing
    for i in `seq 1 22`; do
        if [ ! -f "${BEDNAME}_filtered-updated-chr${i}.check.log" ]; then
            python ${checkVCF_HOME}/checkVCF.py -r $genome_LOC -o ${BEDNAME}_filtered-updated-chr${i} ${BEDNAME}_filtered-updated-chr${i}.vcf

            if ! grep -q "No error found by checkVCF.py" ${BEDNAME}_filtered-updated-chr${i}.check.log; then
                echo "Error found in chromosome ${i}. Stopping the loop."
                break
            fi
        else
            echo "VCF file for chromosome ${i} already processed. Skipping."
        fi
    done    
else
    echo "All VCF files are already processed and ready. Skipping VCF processing."
fi

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

# Rename chromosomes for chromosomes 21 and 22, sort, and clean up
for i in `seq 1 22`; do
        if [ ! -f "../${BEDNAME}_filtered-updated-chr${i}.vcf.gz" ]; then
            $PLINK_HOME/bcftools annotate --rename-chrs chr_rename.txt ${BEDNAME}_filtered-updated-chr${i}.vcf -Ov -o ${BEDNAME}_filtered-updated-chr${i}_renamed.vcf
            $PLINK_HOME/vcf-sort ${BEDNAME}_filtered-updated-chr${i}_renamed.vcf | bgzip -c >  ../${BEDNAME}_filtered-updated-chr${i}.vcf.gz
            rm ${BEDNAME}_filtered-updated-chr${i}_renamed.vcf
        else
            echo "VCF file for chromosome ${i} already processed. Skipping."
        fi
done    


# Submit imputation job and cleanup only if not already done
cd ../../
bash ./submit_imputation_job.sh ${BEDNAME}_filtered-updated-chr $API_key $STUDY 1:22

# Cleanup
#rm -rf tmp

