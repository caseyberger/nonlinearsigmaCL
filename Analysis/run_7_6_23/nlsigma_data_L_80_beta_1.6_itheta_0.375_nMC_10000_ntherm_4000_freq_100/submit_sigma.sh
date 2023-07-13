#!/bin/bash
#
##SBATCH --job-name=nlsigma_prelim_tests			# Job name
#SBATCH --mail-type=ALL				# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cberger@smith.edu			# Where to send mail
#SBATCH --partition=phyq			# Which partition to use
#SBATCH --nodes=1			# Number of nodes
#SBATCH --cpus-per-task=28			# Number of threads per task (OpenMP)
#SBATCH --mem=1gb			# Job memory request
##SBATCH --time=05:00:00			# Time limit hrs:min:sec
#SBATCH --output=nlsigma_prelim_tests%j.log			# Standard output
#SBATCH --error=err_nlsigma_prelim_tests%j.log			# Standard error log

pwd; hostname; date

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

/usr/bin/time -v ./nonlinearsigma inputs.txt

date
