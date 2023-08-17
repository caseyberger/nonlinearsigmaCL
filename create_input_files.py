#! usr/bin/env python

import os, string

'''
In this script, you enter the parameters you want in the simulation, as well as 
the name of the executable and the job name.
The script then loops through the parameters you've given and generates the 
appropriate input file, then generates the submit script to run the code. 
Finally, it copies the executable to the appropriate directory.

It's up to you to run the scripts once they're generated.

'''
#beta = 1/g = 1.6
beta = 1.6
#number of steps in thermalization
ntherm = 4000
#number of monte carlo steps
nMC = 10000
#number of steps between samples
freq = 100
#list of values for lattice length L
L_list = [10,40,80,120,180]
#list of values for itheta (as fractions of pi)
itheta_list = [0.0,0.0625,0.125,1875,0.25,0.3125,0.375,0.4375,0.5,0.5636,0.625,0.6875,0.75,0.8125,0.875,0.9375,1.,1.0625,1.125]

script_name = "nonlinearsigma"
job_name = "nlsigma_prelim_tests"
email = "cberger@smith.edu"
num_cpus = 28
partition = "phyq"


def generate_input_file(length,beta,itheta,nMC,ntherm,freq):
	#generate file extension 
	file_ext  = "L_"+str(length)+"_beta_"+str(beta)+"_itheta_"+str(itheta)+"_nMC_"+str(nMC)
	file_ext += "_ntherm_"+str(ntherm)+"_freq_"+str(freq)
	#create directory
	currdir = os.getcwd()
	work_dir = currdir+'/nlsigma_data_'+file_ext
	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
	'''
		WRITE THE INPUT FILE
	'''
	filename = work_dir+"/inputs.txt"
	input_file = open(filename,'w')
	input_file.write("L = "+str(length)+'\n')
	input_file.write("beta = "+str(beta)+'\n')
	input_file.write("itheta = "+str(itheta)+'\n')
	input_file.write("ntherm = "+str(ntherm)+'\n')
	input_file.write("nMC = "+str(nMC)+'\n')
	input_file.write("freq = "+str(freq)+'\n')
	input_file.close()
	return file_ext

def copy_executable(script_name, file_ext):
	currdir = os.getcwd()
	work_dir = currdir+'/nlsigma_data_'+file_ext
	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
	source_path = currdir+"/"+script_name
	destination_path = work_dir+"/"
	copy_cmd = "cp "+source_path+" "+destination_path
	os.system(copy_cmd)

def generate_slurm_script(script_name,file_ext,job_name,email, partition, cpus_per_task=38):
	#generates the sbatch file that you run with sbatch filename
	#should be paired with the appropriate input file somehow...
	filename = "submit_sigma.sh"
	#create directory (if it doesn't already exist, which it should)
	currdir = os.getcwd()
	work_dir = currdir+'/nlsigma_data_'+file_ext
	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
	slurm_file = open(work_dir+'/'+filename,'w')
	slurm_file.write("#!/bin/bash\n#\n#")
	slurm_file.write("#SBATCH --job-name="+str(job_name)+"\t\t\t# Job name\n")
	slurm_file.write("#SBATCH --mail-type=ALL\t\t\t\t# Mail events (NONE, BEGIN, END, FAIL, ALL)\n")
	slurm_file.write("#SBATCH --mail-user="+str(email)+"\t\t\t# Where to send mail\n")
	slurm_file.write("#SBATCH --partition="+str(partition)+"\t\t\t# Which partition to use\n")
	slurm_file.write("#SBATCH --nodes=1\t\t\t# Number of nodes\n")
	slurm_file.write("#SBATCH --cpus-per-task="+str(cpus_per_task)+"\t\t\t# Number of threads per task (OpenMP)\n")
	slurm_file.write("#SBATCH --mem=1gb\t\t\t# Job memory request\n")
	slurm_file.write("##SBATCH --time=05:00:00\t\t\t# Time limit hrs:min:sec\n")
	slurm_file.write("#SBATCH --output="+str(job_name)+"%j.log\t\t\t# Standard output\n")
	slurm_file.write("#SBATCH --error=err_"+str(job_name)+"%j.log\t\t\t# Standard error log\n\n")
	slurm_file.write("pwd; hostname; date\n\n")
	slurm_file.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n\n")
	slurm_file.write("/usr/bin/time -v ./nonlinearsigma inputs.txt\n\n")
	slurm_file.write("date\n")
	slurm_file.close()



for length in L_list:
	for itheta in itheta_list:
		file_ext = generate_input_file(length,beta,itheta,nMC,ntherm,freq)
		generate_slurm_script(script_name,file_ext,job_name,email,partition,cpus_per_task = num_cpus)
		copy_executable(script_name, file_ext)