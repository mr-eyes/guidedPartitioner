#!/usr/bin/env bash
#! /bin/bash -login
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mabuelanin@gmail.com
#SBATCH -p bmh
#SBATCH -J genesPartitions_assembly_plass
#SBATCH --time=30:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=50gb
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err


# activate conda in general (MUST BE PLACE BEFORE ANYTHING)
. /home/mhussien/miniconda3/etc/profile.d/conda.sh
conda activate omnigraph


# exit when any command fails
set -e
set -o nounset
set -o errexit
set -x
set -euox pipefail


# make a directory specific to user and job
export MYTMP=/scratch/${USER}/slurm_${SLURM_JOB_NAME}_${SLURM_JOBID}
mkdir -p "$MYTMP"

# force clean it up
function cleanup() {
    cp -r $MYTMP "$SLURM_SUBMIT_DIR"
    rm -rf "$MYTMP";
}

trap cleanup EXIT

# ~~~~~~~~~~~~~~~~~~~~~ START Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd "$MYTMP"

## Copy gene partitions
cp -r /home/mhussien/guidedPartitioner/results/genes_partitions ./genes_partitions

# Interlaced fasta extraction script
cp /home/mhussien/guidedPartitioner/extract_interlaced.py ./


# Download PLASS linux static version
wget https://mmseqs.com/plass/plass-static_sse41.tar.gz; tar xvfz plass-static_sse41.tar.gz
PLASS=./plass/bin/plass

# Set Global Variables
THREADS=32
TMP_DIR=tmp_plass
GENES_PARTITIONS=genes_partitions


############################## (1) START Assembly ####################################

OUTPUT=assembled_transcripts
mkdir ${OUTPUT}

TOTAL_FILES_NUMBER=$(ls -1q ${OUTPUT}/* | wc -l)

echo -e "Processing ${TOTAL_FILES_NUMBER} Partitions...\n"

PLASS_LOG=plass_assembly.log
MERGED_TRANSCRIPTS=all_transcripts.fa

touch ${MERGED_TRANSCRIPTS}

for FASTA in ${GENES_PARTITIONS}/*
do

    # Extract the interlaced Fasta and convert to Fastq
    GENE_ID=$(echo "$FASTA" | cut -f 1 -d '.')
    python extract_interlaced.py ${FASTA}
    R1=${GENE_ID}_1.fastq
    R2=${GENE_ID}_2.fastq

    ${PLASS} nuclassemble -v 0 ${R1} ${R2} assembled_${GENE_ID}.fa ${TMP_DIR} --threads ${THREADS} &> ${PLASS_LOG}

    # Move the output assembled transcripts into ${OUTPUT}
    cat assembled_${GENE_ID}.fa >> ${MERGED_TRANSCRIPTS}
    mv assembled_${GENE_ID}.fa ${OUTPUT}
    
    # Remove temporary files
    rm -rf ${TMP_DIR}
    rm -rf ${R1} ${R2}

    # Print remaining partitions
    TOTAL_FILES_NUMBER=`expr $TOTAL_FILES_NUMBER - 1`
    echo "Remaining partitions: ${TOTAL_FILES_NUMBER}"

done




############################## DONE Assembly #######################################

rm -rf ${GENES_PARTITIONS}

########################### Move scratch to home dir ###################


cp -r $MYTMP "$SLURM_SUBMIT_DIR"
scontrol show job $SLURM_JOB_ID     # Print out final statistics about resource uses before job exits

# Print out values of the current jobs SLURM environment variables
env | grep SLURM

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch