#!/usr/bin/env bash
#! /bin/bash -login
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mabuelanin@gmail.com
#SBATCH -p bmh
#SBATCH -J wholeData_assembly_plass
#SBATCH --time=30:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16gb
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

## Copy raw data
cp /home/mhussien/omnigraph/data/SRR11015356_1.fasta ./
cp /home/mhussien/omnigraph/data/SRR11015356_2.fasta ./


# Convert fasta to fastq
cp /home/mhussien/guidedPartitioner/fasta_to_fastq.py ./

python fasta_to_fastq.py SRR11015356_1.fasta SRR11015356_2.fasta

rm -rf SRR11015356_1.fasta
rm -rf SRR11015356_2.fasta

R1=SRR11015356_1.fastq
R2=SRR11015356_2.fastq


# Set Global Variables
THREADS=4
OUTPUT_DIR=tmp_SRR11015356_plass


############################## (1) START Dumping ####################################


/usr/bin/time -v plass nuclassemble -v 3 ${R1} ${R2} assembled_SRR11015356.fa ${OUTPUT_DIR} --threads ${THREADS}


############################## DONE Partitioning #######################################


########################### Move scratch to home dir ###################


cp -r $MYTMP "$SLURM_SUBMIT_DIR"
scontrol show job $SLURM_JOB_ID     # Print out final statistics about resource uses before job exits

# Print out values of the current jobs SLURM environment variables
env | grep SLURM

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch