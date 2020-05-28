#!/usr/bin/env bash
#! /bin/bash -login
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mabuelanin@gmail.com
#SBATCH -p bmh
#SBATCH -J assessment_genePartitions_rnaQuast
#SBATCH --time=10:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=100gb
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
  rm -rf "$MYTMP"
}

trap cleanup EXIT

# ~~~~~~~~~~~~~~~~~~~~~ START Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd "$MYTMP"

## Copy reference genome db
cp /home/mhussien/guidedPartitioner/data/mouse_genome/* ./

## Copy plass assembled transcripts
cp /home/mhussien/guidedPartitioner/slurm_jobs/genesPartitionsAssembly/minContig1000/postAssembly/slurm_postAssembly_reprsExtraction_21980554/clusters_0.97_cdhit_097_geneParts_assembled ./assembled_transcripts.fa

# Set Global Variables
THREADS=32
GTF=gencode.vM25.chr_patch_hapl_scaff.annotation.gtf
REF_GENOME=GRCm38.p6.genome.fa
TRANSCRIPTS=assembled_transcripts.fa

############################## (1) START Dumping ####################################

/usr/bin/time -v rnaQUAST.py --disable_infer_transcripts --disable_infer_genes --transcripts ${TRANSCRIPTS} --reference ${REF_GENOME} --gtf ${GTF} -t ${THREADS} -o rnaQuast_plass_SRR11015356

############################## DONE Partitioning #######################################

########################### Move scratch to home dir ###################

cp -r $MYTMP "$SLURM_SUBMIT_DIR"
scontrol show job $SLURM_JOB_ID # Print out final statistics about resource uses before job exits

# Print out values of the current jobs SLURM environment variables
env | grep SLURM

scontrol show job ${SLURM_JOB_ID} # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
