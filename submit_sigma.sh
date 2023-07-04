#!/bin/bash
#SBATCH --job-name=nonlinearsigma_test           # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cberger@smith.edu            # Where to send mail	
#SBATCH --partition=phyq                         # Which partition to use
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --mem=1gb                                # Job memory request
##SBATCH --time=05:00:00                          # Time limit hrs:min:sec
#SBATCH --output=nonlinearsigma_test_%j.log      # Standard output 
#SBATCH --error=err_nonlinearsigma_test_%j.log   # Standard output and error log

pwd; hostname; date

echo "Running nonlinear sigma on single CPU core"

./nonlinearsigma inputs.txt

date
