#!/usr/bin/env bash
#! /bin/bash -login
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mabuelanin@gmail.com
#SBATCH -p bmh
#SBATCH -J postAssembly_reprsExtraction
#SBATCH --time=10:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=150gb
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

## Copy assembled genes partitions
cp /home/mhussien/guidedPartitioner/slurm_jobs/genesPartitionsAssembly/minContig1000/slurm_genesPartitions_assembly_plass_21902325/all_transcripts.fa ./all_transcripts.fa


# Set Variables
THREADS=32
MAX_RAM_MB=120000

############################## (1) START Assembly ####################################

OUT_PREFIX=cdhit_097_geneParts_assembled

WORD_SIZE=11
SIM=0.97
REF_FASTA=all_transcripts.fa

cd-hit-est -i ${REF_FASTA} -n ${WORD_SIZE} -c ${SIM} -o clusters_${SIM}_${OUT_PREFIX} -d 0 -T ${THREADS} -M ${MAX_RAM_MB}

echo "Extracting representatives from clusters at %${SIM}"
cat clusters_${SIM}_${OUT_PREFIX}.clstr | grep "\*" | awk -F"[>.]" '{print ">"$2}' | grep -Fwf - -A1 <(seqkit seq -w 0 ${REF_FASTA}) | grep -v "^\-\-" > reps_unitigs_${OUT_PREFIX}_${SIM}.fa

############################## DONE Assembly #######################################

########################### Move scratch to home dir ###################


cp -r $MYTMP "$SLURM_SUBMIT_DIR"
scontrol show job $SLURM_JOB_ID     # Print out final statistics about resource uses before job exits

# Print out values of the current jobs SLURM environment variables
env | grep SLURM

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch