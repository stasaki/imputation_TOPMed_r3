conda activate gwas_imputation

# List of directories to process
directories=("San1_CHOP_RMM_AA" "BU_genotyping" "San1_Broad_1709" "San1_CHOP_382" "Diverse_Joint" "AMPAD_WGS" "ADGC_RMM_AA")

# Create an array to store the names of individual chromosome files
chrom_files=()
chnum=1
TMPDIR="./tmp_duplicate/"
mkdir $TMPDIR
for dir in "${directories[@]}"; do

  # Add the resulting binary files to the list
  chrom_files+=("$TMPDIR/${dir}-chr$chnum")
  continue
  
  path_to_vcf=$(find "$dir" -type f -path "*/job-*/local/chr${chnum}.dose.vcf.gz")
  
  # Convert VCF to PLINK binary format
  plink --vcf ${path_to_vcf} --keep-allele-order --make-bed --out $TMPDIR/s1_chr$chnum
  
  # Perform a self-merge to identify problematic SNPs
  plink --bfile $TMPDIR/s1_chr$chnum --bmerge $TMPDIR/s1_chr$chnum --merge-mode 6 --out $TMPDIR/plink_chr$chnum
  
  # Exclude SNPs that were identified as problematic during the self-merge
  plink --bfile $TMPDIR/s1_chr$chnum --exclude $TMPDIR/plink_chr$chnum.missnp --keep-allele-order --make-bed --out $TMPDIR/s2_chr$chnum
  
  # List duplicate variants
  plink --bfile $TMPDIR/s2_chr$chnum --list-duplicate-vars --out $TMPDIR/plink_chr$chnum
  
  # Exclude duplicate variants
  plink --bfile $TMPDIR/s2_chr$chnum --exclude $TMPDIR/plink_chr$chnum.dupvar --keep-allele-order --make-bed --out $TMPDIR/${dir}-chr$chnum
  
  # Rename IID and FID in the .fam file
  awk -v dir="${dir}" '{print dir "-"$1, dir "-"$2, $3, $4, $5, $6}' $TMPDIR/${dir}-chr$chnum.fam > $TMPDIR/${dir}-$chnum.fam.tmp
  mv $TMPDIR/${dir}-$chnum.fam.tmp $TMPDIR/${dir}-chr$chnum.fam
  
  # Clean up intermediate files for the current chromosome
  rm $TMPDIR/s1_chr$chnum.* $TMPDIR/s2_chr$chnum.* $TMPDIR/plink_chr$chnum.*
done

# Merge all the resulting binary files
merge_list_file=$TMPDIR/merge_list.txt
printf "%s\n" "${chrom_files[@]:1}" > $merge_list_file  # Exclude the first file from the merge list

# Base file for merging
base_file=${chrom_files[0]}

plink --bfile $base_file --merge-list $merge_list_file --keep-allele-order --make-bed --out $TMPDIR/merged_data


# Exclude problematic SNPs from the dataset
plink --bfile $base_file --exclude $TMPDIR/merged_data-merge.missnp --keep-allele-order --make-bed --out $TMPDIR/base_file_clean
for file in "${chrom_files[@]:1}"; do
  plink --bfile $file --exclude $TMPDIR/merged_data-merge.missnp --keep-allele-order --make-bed --out $TMPDIR/$(basename $file)_clean
done

# Create a new merge list with cleaned files
clean_files=()
for file in "${chrom_files[@]}"; do
  clean_files+=("$TMPDIR/$(basename $file)_clean")
done
printf "%s\n" "${clean_files[@]:1}" > $merge_list_file  # Exclude the first file from the merge list

# Merge the cleaned files
plink --bfile $TMPDIR/base_file_clean --merge-list $merge_list_file --keep-allele-order --make-bed --out $TMPDIR/merged_data_clean


# Filter SNPs with too many missing values
plink --bfile $TMPDIR/merged_data_clean --geno 0.05 --maf 0.05 --make-bed --out $TMPDIR/filtered_data

# Prune the filtered data
plink --bfile $TMPDIR/filtered_data --indep-pairwise 50 5 0.2 --out $TMPDIR/pruned_data

# Compute IBD on pruned data
plink --bfile $TMPDIR/filtered_data --extract $TMPDIR/pruned_data.prune.in --genome --out $TMPDIR/ibd_results

# Optionally, remove the merge list file if no longer needed
rm $merge_list_file

# run duplicte_check.R

rm -r $TMPDIR