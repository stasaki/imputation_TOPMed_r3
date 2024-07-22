#bash 01_get_unique_projids.sh 

# https://radc.rush.edu/radc/dynamic/reports/projidStudySummaryByList.htm
# validate projids
# saved as do_not_distribute.txt

wget https://www.chg.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
#The TOPMed reference panel is not available for direct download from this site, it needs to be created from the VCF of dbSNP submitted sites (currently ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz).
#This can be downloaded from the Bravo Website https://bravo.sph.umich.edu/freeze5/hg38/
wget https://www.chg.ox.ac.uk/~wrayner/tools/CreateTOPMed.zip
./CreateTOPMed.pl -i bravo-dbsnp-all.vcf.gz

bash imputation.sh gs://shinya_test/Genotype/AMPAD_WGS/03_liftover/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08.recalibrated_variants_hg38_dn.QCed1 AMPAD_WGS
bash imputation.sh gs://shinya_test/RADC_omics/WGS/03_Diverse_Joint/DIVERSE.joint.callset Diverse_Joint


cd AMPAD_WGS
nohup /data/shinya/shinya/Resource/GENETICS_Resource/softwares/imputationbot/imputationbot  download  job-20240711-000420-579 --password ""

cd Diverse_Joint
nohup /data/shinya/shinya/Resource/GENETICS_Resource/softwares/imputationbot/imputationbot  download  job-20240711-011328-645 --password ""

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



# Set temporary directory
TMPDIR=tmp
mkdir -p $TMPDIR

# Loop through chromosomes
for chr in {22..1}; do
    merged_vcf="$TMPDIR/merged_chr${chr}.vcf.gz"

    # Check if the merged VCF file already exists
    if [ -f "$merged_vcf" ]; then
        echo "Merged VCF for chromosome $chr already exists. Skipping."
        continue
    fi

    echo "Processing chromosome $chr..."

    # Filter and annotate Diverse VCF
    bcftools view -S ^$TMPDIR/Diverse_remove.txt -O u "Diverse_Joint/job-20240711-011328-645-Diverse_Joint/local/chr${chr}.dose.vcf.gz" | \
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O z -o $TMPDIR/div.vcf.gz
    bcftools index $TMPDIR/div.vcf.gz
    
    # Annotate WGS VCF
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O z -o $TMPDIR/wgs.vcf.gz AMPAD_WGS/job-20240711-000420-579-AMPAD_WGS/local/chr${chr}.dose.vcf.gz
    bcftools index $TMPDIR/wgs.vcf.gz
    
    # Filter and annotate Array VCF
    bcftools view -S ^$TMPDIR/Array_remove.txt -O u "combined_vcf/merged_chr${chr}.vcf.gz" | \
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O z -o $TMPDIR/array.vcf.gz
    bcftools index $TMPDIR/array.vcf.gz
    
    # Convert VCF files to PLINK format
    plink --vcf $TMPDIR/wgs.vcf.gz --keep-allele-order --make-bed --out $TMPDIR/wgs
    plink --vcf $TMPDIR/array.vcf.gz --keep-allele-order --make-bed --out $TMPDIR/array
    plink --vcf $TMPDIR/div.vcf.gz --keep-allele-order --make-bed --out $TMPDIR/div

    # Create a merge list for PLINK
    merge_list="$TMPDIR/merge_list.txt"
    echo "$TMPDIR/array" > "$merge_list"
    echo "$TMPDIR/div" >> "$merge_list"

    # Merge PLINK files
    plink --bfile $TMPDIR/wgs --merge-list $merge_list --geno 0.1 --keep-allele-order --make-bed --out $TMPDIR/combined

    # Extract IDs from merged PLINK file
    cut -f2 $TMPDIR/combined.bim > $TMPDIR/keep.id

    # Create a list of VCF files for merging
    vcf_list_file="$TMPDIR/chr${chr}_vcf_list.txt"
    echo "$TMPDIR/div.vcf.gz" > "$vcf_list_file"
    echo "$TMPDIR/array.vcf.gz" >> "$vcf_list_file"
    echo "$TMPDIR/wgs.vcf.gz" >> "$vcf_list_file"

    # Check if any VCF files were found
    if [ ! -s "$vcf_list_file" ]; then
        echo "No VCF files found for chromosome $chr. Skipping."
        continue
    fi

    # Filter each VCF file and save the filtered versions
    while read vcf; do
        filtered_vcf="${vcf%.vcf.gz}_filtered.vcf.gz"
        bcftools view -i "ID=@$TMPDIR/keep.id" -O z -o "$filtered_vcf" "$vcf"
        bcftools index $filtered_vcf
    done < "$vcf_list_file"

    # Create a list of filtered VCF files
    filtered_vcf_list_file="$TMPDIR/filtered_vcf_list.txt"
    ls $TMPDIR/*_filtered.vcf.gz > "$filtered_vcf_list_file"

    # Merge the filtered VCF files
    bcftools merge -m both -O z -l "$filtered_vcf_list_file" -o "$merged_vcf"

    # Clean up temporary files
    rm "$vcf_list_file" $TMPDIR/div.vcf.gz $TMPDIR/div.vcf.gz.csi $TMPDIR/array.vcf.gz $TMPDIR/array.vcf.gz.csi
done

VCF="./tmp"
TMPDIR="./tmp_plink"

# Create an array to store the names of individual chromosome files
chrom_files=()

for chnum in {1..22}; do
  # Convert VCF to PLINK binary format
  plink --vcf $VCF/merged_chr$chnum.vcf.gz --keep-allele-order --make-bed --out $TMPDIR/s1_chr$chnum
  
  ## Perform a self-merge to identify problematic SNPs
  #plink --bfile $TMPDIR/s1_chr$chnum --bmerge $TMPDIR/s1_chr$chnum --merge-mode 6 --out $TMPDIR/plink_chr$chnum
  
  ## Exclude SNPs that were identified as problematic during the self-merge
  #plink --bfile $TMPDIR/s1_chr$chnum --exclude $TMPDIR/plink_chr$chnum.missnp --keep-allele-order --make-bed --out $TMPDIR/s2_chr$chnum
  
  ## List duplicate variants
  #plink --bfile $TMPDIR/s2_chr$chnum --list-duplicate-vars --out $TMPDIR/plink_chr$chnum
  
  ## Exclude duplicate variants
  #plink --bfile $TMPDIR/s2_chr$chnum --exclude $TMPDIR/plink_chr$chnum.dupvar --keep-allele-order --make-bed --out $TMPDIR/s3_chr$chnum

  # Add the final chromosome file to the array
  chrom_files+=($TMPDIR/s1_chr$chnum)
  
  # Clean up intermediate files for the current chromosome
  #rm $TMPDIR/s1_chr$chnum.* $TMPDIR/s2_chr$chnum.* $TMPDIR/plink_chr$chnum.*
done

# Merge all chromosomes into a single PLINK binary file
# Create a temporary file to hold the merge list
merge_list=$TMPDIR/merge_list.txt
for ((i = 1; i < ${#chrom_files[@]}; i++)); do
  echo ${chrom_files[i]} >> $merge_list
done

plink --bfile ${chrom_files[0]} --merge-list $merge_list --keep-allele-order --make-bed --out $TMPDIR/combined

#/home/shinya/miniconda3/envs/reva_5hmc/bin/python rename_variants_with_backup.py ./tmp_plink/*.bim --problematic ./tmp_plink/combined-merge.missnp

mv tmp_plink/combined.bim combined_plink/TOPMed_r3_all_hg38.bim
mv tmp_plink/combined.bed combined_plink/TOPMed_r3_all_hg38.bed
mv tmp_plink/combined.fam combined_plink/TOPMed_r3_all_hg38.fam



plink --bfile TOPMed_r3_all_hg38 --recode vcf --out intermediate_output

awk '{print $1 "_" $2, $2}'  TOPMed_r3_all_hg38.fam   > sample_map.txt
 
 # Compress the VCF file
bgzip -c intermediate_output.vcf >intermediate_output.vcf.gz


# Rename samples in the compressed VCF
bcftools reheader -s sample_map.txt -o TOPMed_r3_all_hg38.vcf.gz intermediate_output.vcf.gz

# Index the compressed VCF file
bcftools index TOPMed_r3_all_hg38.vcf.gz
tabix  TOPMed_r3_all_hg38.vcf.gz

rm sample_map.txt 
rm intermediate_output.*



