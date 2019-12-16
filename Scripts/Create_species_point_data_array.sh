#!/bin/bash

#SBATCH --partition=milkun
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=sp_select
#SBATCH --mem-per-cpu=32G
#SBATCH --array=1-1936%100
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o "/vol/milkunB/Projects/Bioclimatic_envelopes/Output/Logs/Processing/log_create_pt_data_%a.out"
#SBATCH --mail-type=FAIL


export TMPDIR=/scratch/bclim
mkdir -p $TMPDIR

srun /opt/R-3.4.2/bin/R --vanilla --no-save --args ${SLURM_ARRAY_TASK_ID} < /vol/milkunB/Projects/Bioclimatic_envelopes/Analyses/Processing/Create_species_point_data.R
