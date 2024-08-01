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
bash ./script/imputation.sh gs://shinya_test/Genotype/ADGC_RMM_AA/adgc.rosmapmars.afam.vFinal_uniq_pos_hg38 ADGC_RMM_AA

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

find . -type f -path "./combined_vcf/*.vcf.gz" | while read vcf_file; do
    index_file="${vcf_file}.csi"
    if [ ! -f "$index_file" ]; then
        echo "Indexing: $vcf_file"
        bcftools index "$vcf_file"
    else
        echo "Index already exists for: $vcf_file"
    fi
done



VCF="./tmp"
TMPDIR="./tmp_plink"

mkdir $VCF
mkdir $TMPDIR

# Function to process a single chromosome
process_chromosome() {
    local chr=$1
    echo "Merging chromosome $chr..."

    # Create a temporary file to store VCF paths
    vcf_list_file="./tmp/chr${chr}_vcf_list.txt"

    # Find and list all VCF files for the current chromosome
    #for dir in "${directories[@]}"; do
    #echo $dir
    #    find "$dir" -type f -path "*/job-*/local/chr${chr}.dose.vcf.gz" >> "$vcf_list_file"
    #done
    echo "./ADGC_RMM_AA/job-20240725-015404-612-ADGC_RMM_AA/local/chr${chr}.dose.vcf.gz" >> "$vcf_list_file"
    echo "./combined_vcf/TOPMed_r3_array_chr${chr}.vcf.gz" >> "$vcf_list_file"

    # Check if any VCF files were found
    if [ ! -s "$vcf_list_file" ]; then
        echo "No VCF files found for chromosome $chr. Skipping."
        return
    fi

    # Merge VCF files
    merged_vcf="./tmp/merged_chr${chr}.vcf.gz"
    bcftools merge -m both -O z -o "$merged_vcf" -l "$vcf_list_file"

    # Remove temporary file
    rm "$vcf_list_file"
}

export -f process_chromosome

# Run the function in parallel for chromosomes 1 to 22
seq 1 22 | xargs -n 1 -P 3 bash -c 'process_chromosome "$@"' _ 


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


plink --bfile ${chrom_files[0]} --merge-list $merge_list --keep-allele-order --make-bed --out $TMPDIR/combined
/home/shinya/miniconda3/envs/reva_5hmc/bin/python ./script/rename_variants_with_backup.py ./tmp_plink/*.bim --problematic ./tmp_plink/combined-merge.missnp
plink --bfile ${chrom_files[0]} --merge-list $merge_list --keep-allele-order --make-bed --out $TMPDIR/combined

mv $TMPDIR/combined.bim combined_plink/TOPMed_r3_array_hg38.bim
mv $TMPDIR/combined.bed combined_plink/TOPMed_r3_array_hg38.bed
mv $TMPDIR/combined.fam combined_plink/TOPMed_r3_array_hg38.fam

for chnum in {1..22}; do
  mv $VCF/merged_chr${chnum}.vcf.gz ./combined_vcf/TOPMed_r3_array_chr${chnum}.vcf.gz
done
rm -r $VCF
rm -r $TMPDIR

###############
# Set temporary directory
TMPDIR=tmp
mkdir -p $TMPDIR

# Function to process a single chromosome
process_chromosome() {
    local chr=$1
    TMPDIR=tmp2
    merged_vcf="$TMPDIR/merged_chr${chr}.vcf.gz"

    # Check if the merged VCF file already exists
    if [ -f "$merged_vcf" ]; then
        echo "Merged VCF for chromosome $chr already exists. Skipping."
        return
    fi

    echo "Processing chromosome $chr..."

    # Filter and annotate Diverse VCF
    bcftools view -S ^$TMPDIR/ADGC_RMM_AA_remove.txt -O u "ADGC_RMM_AA/job-20240725-015404-612-ADGC_RMM_AA/local/chr${chr}.dose.vcf.gz" | \
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O z -o $TMPDIR/ADGC_${chr}.vcf.gz
    bcftools index $TMPDIR/ADGC_${chr}.vcf.gz
    
    # Annotate all VCF
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O z -o $TMPDIR/all_${chr}.vcf.gz combined_vcf/TOPMed_r3_all_chr${chr}.vcf.gz
    bcftools index $TMPDIR/all_${chr}.vcf.gz

    # Convert VCF files to PLINK format
    #plink --vcf $TMPDIR/ADGC_${chr}.vcf.gz --keep-allele-order --make-bed --out $TMPDIR/ADGC_${chr}
    plink --vcf $TMPDIR/all_${chr}.vcf.gz --keep-allele-order --make-bed --out $TMPDIR/all_${chr}

    # Create a merge list for PLINK
    #merge_list="$TMPDIR/merge_list_${chr}.txt"
    #echo "$TMPDIR/ADGC_${chr}" > "$merge_list"
    #echo "$TMPDIR/all_${chr}" >> "$merge_list"

    # Merge PLINK files
    #plink --bfile $TMPDIR/all_${chr} --merge-list $merge_list --geno 0.1 --keep-allele-order --make-bed --out $TMPDIR/combined_all_${chr}

    # Extract IDs from merged PLINK file
    cut -f2 $TMPDIR/all_${chr}.bim > $TMPDIR/keep.id_${chr}

    bcftools view -i "ID=@$TMPDIR/keep.id_${chr}" -O z -o "$TMPDIR/ADGC_${chr}_filtered.vcf.gz" "$TMPDIR/ADGC_${chr}.vcf.gz"
    bcftools index $TMPDIR/ADGC_${chr}_filtered.vcf.gz

    # Create a list of filtered VCF files
    filtered_vcf_list_file="$TMPDIR/chr${chr}_filtered_vcf_list.txt"
    echo $TMPDIR/ADGC_${chr}_filtered.vcf.gz > "$filtered_vcf_list_file"
    echo $TMPDIR/all_${chr}.vcf.gz >> "$filtered_vcf_list_file"
    
    # Merge the filtered VCF files
    bcftools merge -m both -O z -l "$filtered_vcf_list_file" -o "$merged_vcf"

    # Clean up temporary files
    rm $TMPDIR/ADGC_${chr}.vcf.gz $TMPDIR/ADGC_${chr}.vcf.gz.csi $TMPDIR/all_${chr}.vcf.gz $TMPDIR/all_${chr}.vcf.gz.csi
    rm $TMPDIR/all_${chr}.*
    rm $TMPDIR/keep.id_${chr}
    rm $TMPDIR/ADGC_${chr}_filtered.vcf.gz $TMPDIR/ADGC_${chr}_filtered.vcf.gz.csi 
    rm "$TMPDIR/chr${chr}_filtered_vcf_list.txt"
}

export -f process_chromosome

# Run the function in parallel for chromosomes 1 to 22
seq 1 22 | xargs -n 1 -P 3 bash -c 'process_chromosome "$@"' _ 


VCF="./tmp"
TMPDIR="./tmp_plink"

# Create an array to store the names of individual chromosome files
chrom_files=()

for chnum in {1..22}; do
  # Convert VCF to PLINK binary format
  plink --vcf $VCF/merged_chr$chnum.vcf.gz --keep-allele-order --make-bed --out $TMPDIR/s1_chr$chnum

  # Add the final chromosome file to the array
  chrom_files+=($TMPDIR/s1_chr$chnum)
  
done

# Merge all chromosomes into a single PLINK binary file
# Create a temporary file to hold the merge list
merge_list=$TMPDIR/merge_list.txt
for ((i = 1; i < ${#chrom_files[@]}; i++)); do
  echo ${chrom_files[i]} >> $merge_list
done


plink --bfile ${chrom_files[0]} --merge-list $merge_list --keep-allele-order --make-bed --out $TMPDIR/combined

mv $TMPDIR/combined.bim combined_plink/TOPMed_r3_all_hg38.bim
mv $TMPDIR/combined.bed combined_plink/TOPMed_r3_all_hg38.bed
mv $TMPDIR/combined.fam combined_plink/TOPMed_r3_all_hg38.fam

cd combined_plink
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

cd ..
for chnum in {1..22}; do
  mv $VCF/merged_chr${chnum}.vcf.gz ./combined_vcf/TOPMed_r3_all_chr${chnum}.vcf.gz
done
rm -r $VCF
rm -r $TMPDIR