#!/bin/bash

# Function to process chromosome range or list
process_chromosomes() {
    local RANGE_REGEX='^([0-9]+):([0-9]+)$'
    local CHR_NUMS=()

    for ARG in "$@"; do
        if [[ $ARG =~ $RANGE_REGEX ]]; then
            for ((i=${BASH_REMATCH[1]}; i<=${BASH_REMATCH[2]}; i++)); do
                CHR_NUMS+=($i)
            done
        else
            CHR_NUMS+=($ARG)
        fi
    done

    echo "${CHR_NUMS[@]}"
}

# Check if the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <path-to-files> <token> <job-name> <chr_range_or_list>"
  echo "Example: $0 /path-to mytoken my_job_name 1:22"
  echo "Example: $0 /path-to mytoken my_job_name 1 2 21 22"
  exit 1
fi

# Extract the path-to, token, and job-name arguments
PATH_TO_FILES=$1
TOKEN=$2
JOB_NAME=$3

# Process chromosomes (range or list)
shift 3 # Skip the first three arguments
CHROMOSOMES=$(process_chromosomes "$@")

# Construct the files part of the curl command
FILES_PART=""
for CHR_NUM in $CHROMOSOMES; do
    FILES_PART+="-F \"files=@${PATH_TO_FILES}${CHR_NUM}.vcf.gz\" "
done

# Construct the full curl command
CURL_COMMAND="curl https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver \
  -X \"POST\" \
  -H \"X-Auth-Token: ${TOKEN}\" \
  ${FILES_PART}\
  -F \"job-name=${JOB_NAME}\" \
  -F \"refpanel=apps@topmed-r3\" \
  -F \"build=hg38\" \
  -F \"phasing=eagle\" \
  -F \"population=all\" \
  -F \"meta=yes\" \
  -F \"mode=imputation\" \
  -F \"r2Filter=0.3\""

# Execute the curl command
#eval $CURL_COMMAND
echo $CURL_COMMAND
