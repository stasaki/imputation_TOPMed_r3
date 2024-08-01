#bash 01_get_unique_projids.sh 

# https://radc.rush.edu/radc/dynamic/reports/projidStudySummaryByList.htm
# validate projids
# saved as do_not_distribute.txt

wget https://www.chg.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
#The TOPMed reference panel is not available for direct download from this site, it needs to be created from the VCF of dbSNP submitted sites (currently ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz).
#This can be downloaded from the Bravo Website https://bravo.sph.umich.edu/freeze5/hg38/
wget https://www.chg.ox.ac.uk/~wrayner/tools/CreateTOPMed.zip
./CreateTOPMed.pl -i bravo-dbsnp-all.vcf.gz
 git clone https://github.com/vcftools/vcftools.git
 
conda activate gwas_imputation
bash ./script/imputation.sh gs://shinya_test/RADC_omics/Genotype/San1_Broad_1709/03_imputation/output/unimputed_bed/Merge.phaseIV.final_uniq_pos_QCed_hg38 San1_Broad_1709
bash ./script/imputation.sh gs://shinya_test/RADC_omics/Genotype/San1_CHOP_382/03_imputation/output/unimputed_bed/chop.rosmap.euam.vFinal_uniq_pos_QCed_hg38 San1_CHOP_382
bash ./script/imputation.sh gs://shinya_test/RADC_omics/Genotype/BU_genotyping/03_imputation/output/unimputed_bed/BU_measured_uniq_pos_QCed_hg38 BU_genotyping
bash ./script/imputation.sh gs://shinya_test/RADC_omics/Genotype/San1_CHOP_RMM_AA/03_imputation/output/unimputed_bed/chop.rosmap.afam.vFinal_uniq_pos_QCed_hg38 San1_CHOP_RMM_AA
bash ./script/imputation.sh gs://shinya_test/Genotype/ADGC_RMM_AA/adgc.rosmapmars.afam.vFinal_uniq_pos_hg38 ADGC_RMM_AA


cd San1_Broad_1709
/home/shinya/Resource/GENETICS_Resource/softwares/imputationbot/imputationbot  download  job-20240401-213056-501 --password ""

cd San1_CHOP_382
nohup /home/shinya/Resource/GENETICS_Resource/softwares/imputationbot/imputationbot  download  job-20240401-212844-400 --password ""

cd BU_genotyping
nohup /home/shinya/Resource/GENETICS_Resource/softwares/imputationbot/imputationbot  download  job-20240328-012610-140 --password ""

cd San1_CHOP_RMM_AA
nohup /home/shinya/Resource/GENETICS_Resource/softwares/imputationbot/imputationbot  download  job-20240401-203900-049 --password ""


cd ADGC_RMM_AA
nohup /data/shinya/shinya/Resource/GENETICS_Resource/softwares/imputationbot/imputationbot  download  job-20240725-015404-612 --password ""


# Find all VCF files and index them
conda activate gwas_imputation
find . -type f -path "./*/job-*/local/*.dose.vcf.gz" | while read vcf_file; do
    index_file="${vcf_file}.csi"
    if [ ! -f "$index_file" ]; then
        echo "Indexing: $vcf_file"
        bcftools index "$vcf_file"
    else
        echo "Index already exists for: $vcf_file"
    fi
done


# List of directories to process
directories=("San1_CHOP_RMM_AA" "BU_genotyping" "San1_Broad_1709" "San1_CHOP_382")

# Loop through chromosomes
for chr in {1..22}; do
    echo "Merging chromosome $chr..."

    # Create a temporary file to store VCF paths
    vcf_list_file="./tmp/chr${chr}_vcf_list.txt"

    # Find and list all VCF files for the current chromosome
    for dir in "${directories[@]}"; do
        find "$dir" -type f -path "*/job-*/local/chr${chr}.dose.vcf.gz" >> "$vcf_list_file"
    done

    # Check if any VCF files were found
    if [ ! -s "$vcf_list_file" ]; then
        echo "No VCF files found for chromosome $chr. Skipping."
        continue
    fi

    # Merge VCF files
    merged_vcf="./tmp/merged_chr${chr}.vcf.gz"
    bcftools merge -m both -O z -o "$merged_vcf" -l "$vcf_list_file"

    # Remove temporary file
    rm "$vcf_list_file"
done



VCF="./tmp"
TMPDIR="./tmp_plink"

# Create an array to store the names of individual chromosome files
chrom_files=()

for chnum in {1..22}; do
  # Convert VCF to PLINK binary format
  plink --vcf $VCF/merged_chr$chnum.vcf.gz --keep-allele-order --make-bed --out $TMPDIR/s1_chr$chnum
  
  # Perform a self-merge to identify problematic SNPs
  plink --bfile $TMPDIR/s1_chr$chnum --bmerge $TMPDIR/s1_chr$chnum --merge-mode 6 --out $TMPDIR/plink_chr$chnum
  
  # Exclude SNPs that were identified as problematic during the self-merge
  plink --bfile $TMPDIR/s1_chr$chnum --exclude $TMPDIR/plink_chr$chnum.missnp --keep-allele-order --make-bed --out $TMPDIR/s2_chr$chnum
  
  # List duplicate variants
  plink --bfile $TMPDIR/s2_chr$chnum --list-duplicate-vars --out $TMPDIR/plink_chr$chnum
  
  # Exclude duplicate variants
  plink --bfile $TMPDIR/s2_chr$chnum --exclude $TMPDIR/plink_chr$chnum.dupvar --keep-allele-order --make-bed --out $TMPDIR/s3_chr$chnum

  # Add the final chromosome file to the array
  chrom_files+=($TMPDIR/s3_chr$chnum)
  
  # Clean up intermediate files for the current chromosome
  rm $TMPDIR/s1_chr$chnum.* $TMPDIR/s2_chr$chnum.* $TMPDIR/plink_chr$chnum.*
done

# Merge all chromosomes into a single PLINK binary file
# Create a temporary file to hold the merge list
merge_list=$TMPDIR/merge_list.txt
for ((i = 1; i < ${#chrom_files[@]}; i++)); do
  echo ${chrom_files[i]} >> $merge_list
done

/home/shinya/miniconda3/envs/reva_5hmc/bin/python rename_variants_with_backup.py ./tmp_plink/*.bim --problematic ./tmp_plink/combined-merge.missnp


plink --bfile ${chrom_files[0]} --merge-list $merge_list --keep-allele-order --make-bed --out $TMPDIR/combined


# run filter_samples.R
# run extract_and_rename_vars.R

# Define file paths
bedA="tmp_plink/combined"
remove_samples_for_bedA="tmp_plink/Array_remove.txt"
keep_variants_for_bedA="tmp_plink/Array_keep_vars.txt"
bedB="Diverse_Joint/tmp/DIVERSE.joint.callset_filtered-updated"
remove_samples_for_bedB="tmp_plink/Diverse_remove.txt"
keep_variants_for_bedB="tmp_plink/Diverse_keep_vars.txt"
bedC="AMPAD_WGS/tmp/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08.recalibrated_variants_hg38_dn.QCed1_filtered-updated"
keep_variants_for_bedC="tmp_plink/AMPAD_keep_vars.txt"

# Remove specified samples and extract SNPs from bedA
plink --bfile $bedA --remove $remove_samples_for_bedA --extract $keep_variants_for_bedA --make-bed --out tmp_plink/combined_cleaned

# Remove specified samples and extract SNPs from bedB
plink --bfile $bedB --remove $remove_samples_for_bedB --extract $keep_variants_for_bedB --make-bed --out tmp_plink/DIVERSE_cleaned

# Extract SNPs from bedC
plink --bfile $bedC --extract $keep_variants_for_bedC --make-bed --out tmp_plink/AMPAD_cleaned


# Create merge list file
# Create merge list file
echo "tmp_plink/combined_cleaned" > merge_list.txt
echo "tmp_plink/DIVERSE_cleaned" >> merge_list.txt
echo "tmp_plink/AMPAD_cleaned" >> merge_list.txt

# Merge datasets
plink --merge-list merge_list.txt --make-bed --out ./tmp_plink/combined_dataset

plink --bfile ./tmp_plink/combined_dataset --geno 0.1 --make-bed --out ./tmp_plink/filtered_dataset

plink --bfile TOPMed_r3_array_and_WGS_hg38 --recode vcf --out intermediate_output

awk '{print $1 "_" $2, $2}'  TOPMed_r3_array_and_WGS_hg38.fam   > sample_map.txt
 
 # Compress the VCF file
bgzip -c intermediate_output.vcf >intermediate_output.vcf.gz


# Rename samples in the compressed VCF
bcftools reheader -s sample_map.txt -o TOPMed_r3_array_and_WGS_hg38.vcf.gz intermediate_output.vcf.gz

# Index the compressed VCF file
bcftools index TOPMed_r3_array_and_WGS_hg38.vcf.gz

rm TOPMed_r3_array_and_WGS_hg38.vcf
rm TOPMed_r3_array_and_WGS_hg38.log
rm sample_map.txt 
rm intermediate_output.vcf.gz 


