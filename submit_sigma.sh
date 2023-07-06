#!/bin/bash
#SBATCH --job-name=nonlinearsigma_omp_test           # Job name
#SBATCH --mail-type=ALL                              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cberger@smith.edu                # Where to send mail
#SBATCH --partition=phyq                             # Which partition to use
#SBATCH --nodes=1                                    # Number of nodes
#SBATCH --cpus-per-task=28                           # Number of threads per task (OpenMP)
#SBATCH --mem=1gb                                    # Job memory request
##SBATCH --time=05:00:00                             # Time limit hrs:min:sec
#SBATCH --output=nonlinearsigma_omp_test_%j.log      # Standard output 
#SBATCH --error=err_nonlinearsigma_omp_test_%j.log   # Standard output and error log

pwd; hostname; date

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Running nonlinear sigma on single CPU core"

/usr/bin/time -v ./nonlinearsigma inputs.txt

date
