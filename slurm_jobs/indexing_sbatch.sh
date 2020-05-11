#!/usr/bin/env bash
#! /bin/bash -login
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mabuelanin@gmail.com
#SBATCH -p bmh
#SBATCH -J genes_indexing
#SBATCH --time=5:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=60gb
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
    rm -rf "$MYTMP";
}

trap cleanup EXIT

# ~~~~~~~~~~~~~~~~~~~~~ START Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd "$MYTMP"

## Copy raw data
cp /home/mhussien/guidedPartitioner/data/gencode.v34.transcripts.fa ./
cp /home/mhussien/guidedPartitioner/data/gencode.v34.transcripts.fa.names ./


############################## (1) START Indexing ####################################

FASTA=gencode.v34.transcripts.fa
NAMES=gencode.v34.transcripts.fa.names

/usr/bin/time -v /home/mhussien/guidedPartitioner/build/genesIndexing ${FASTA} ${NAMES}

# Remove reads
rm -rf ${FASTA}
rm -rf ${NAMES}

############################## DONE #######################################


########################### Move scratch to home dir ###################


cp -r $MYTMP "$SLURM_SUBMIT_DIR"
scontrol show job $SLURM_JOB_ID     # Print out final statistics about resource uses before job exits

# Print out values of the current jobs SLURM environment variables
env | grep SLURM

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch