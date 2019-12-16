#!/bin/bash

#SBATCH --partition=milkun
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=md_ens
#SBATCH --mem-per-cpu=12G
#SBATCH --array=1-1915%10
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o "/vol/milkunB/Projects/Bioclimatic_envelopes/Output/Logs/Modeling/log_run_ensemble_%a.out"
#SBATCH --mail-type=FAIL

export TMPDIR=/scratch/mdls
mkdir -p $TMPDIR

srun /opt/R-3.4.2/bin/R --vanilla --no-save --args ${SLURM_ARRAY_TASK_ID} < /vol/milkunB/Projects/Bioclimatic_envelopes/Analyses/Modeling/Run_ensemble_models.R
