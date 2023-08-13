# Simulation and Analysis Code for $2D$ $\mathcal{O}(3)$ Nonlinear Sigma Model with Topological Term $\theta$

## Simulation code for Monte Carlo approach

This approach uses a conventional Monte Carlo simulation with a Metropolis step to simulate the $2D$ $\mathcal{O}(3)$ nonlinear sigma model where the topological term $\theta$ is assigned an imaginary value.

### Background and Methodology

#### Lattice Action

In the continuum, the model has action

$S = \frac{1}{2g} \int d^{2} x \left( \partial_{\mu} \vec{\phi}(x)\right)^{2} - i \theta \int d^{2} x Q(x)$

with

$Q(x) = \frac{1}{8 \pi} \epsilon^{\mu \nu} \epsilon_{abc} \partial_{\mu} \phi^{b}(x) \partial_{\nu}\phi^{c}(x)$

with $\phi$ a 3-component unit vector.

In lattice action form, we have (see B. Alles, M. Giordano, and A. Papa. Behavior near θ = π of the mass gap in the two-dimensional o(3) nonlinear sigma model. Phys. Rev. B, 90:184421, Nov 2014):

$
\begin{eqnarray}
S_{L} &=& A_{L} - i \theta Q_{L}\\
A_{L} &=& -\frac{1}{g_{L}}\sum_{x,\mu}\left(\phi_{x}\cdot \phi_{x+\hat{\mu}}\right)\\
Q_{L} &=& \sum_{x}\sum_{\Delta} Q_{L} \Delta
\end{eqnarray}
$

where $\vec{\phi}$ a 3-component unit vector ($\vec{\phi} \cdot\vec{\phi} = 1$) and $Q_{L}$ is the total topological charge on the lattice. 

#### Regularization of the topological charge

The topological charge has been defined via sums over triangles created by cutting each square plaquette along the diagonal. Each vertex is labeled (numbered counter-clockwise), such that we call the fields at the sites of the vertices $\vec{\phi}_{1}$, $\vec{\phi}_{2}$, and $\vec{\phi}_{3}$. 


The topological charge over each triangle obeys

$\exp(2 \pi i Q_{L}(\Delta)) = \frac{1}{\rho}\left(1 + \vec{\phi}_{1}\cdot\vec{\phi}_{2} + \vec{\phi}_{2}\cdot\vec{\phi}_{3} + \vec{\phi}_{3}\cdot\vec{\phi}_{1} + i \vec{\phi}_{1} \cdot (\vec{\phi}_{2}\times\vec{\phi}_{3})\right)$

with 

$\rho^{2} = 2(1+\vec{\phi}_{1}\cdot\vec{\phi}_{2})(1 + \vec{\phi}_{2}\cdot\vec{\phi}_{3})(1+ \vec{\phi}_{3}\cdot\vec{\phi}_{1})$ 

and 

$Q_{L}(\Delta) \in \left[-\frac{1}{2}, \frac{1}{2}\right]$

We use the arcsin of the quantity $\exp(2 \pi i Q_{L}(\Delta))$ to compute $Q_{L}(\Delta)$, as in C++ the domain of arcsin is symmetric about $0$, which prevents the need to adjust the domain to fit the expectation given above.

<div class="alert alert-warning">
When we sum over all unique triangles on the lattice, the total topological charge $Q_{L} = \sum_{\Delta} Q_{}(\Delta)$ should return integer values, however we are currently not finding this to be the case. More investigation is required here.
</div>

#### Analytic continuation

The simulation is run for imaginary values of $\theta$, which means we must analytically continue our results for real $\theta$. To do this, we fit our results to a curve and then substitute our imaginary $i \theta = \nu$ for a real $\theta = -i \nu$.

#### The mass gap

Our goal with this project is to determine the mass gap, which should vanish as $\theta \to \pi$. The mass gap is the inverse of the correlation length, which we compute in the simulation.

<div class="alert alert-warning">
The correlation length should be a real number, but our simulations are currently returning a complex result, due to the complexity of the correlation function. More investigation is required here.
</div>

### Running the simulation code

#### Makefile flags

The code is written in C++ and OpenMP and can be compiled with numerous flags.

```bash
USE_OMP ?= TRUE
USE_GPROF ?= FALSE
USE_TEST_PRINT_STATEMENTS ?= FALSE
USE_EXTREME_TEST_CONDITION ?= FALSE
USE_CHECK_QL_COS ?= FALSE
USE_CONST_RN ?= FALSE
```

The first flag "USE_OMP" toggles whether to implement the parallelization in the code. It should be set to "TRUE" if you want to run the simulation in parallel. This is highly recommended for large lattices as the scaling is very poor in series.

If you wish to profile the code, set the second flag "USE_GPROF" to "TRUE". This sets the correct compiler flags so that you can generate the profiling output. To view the output after the code has run, go to the directory in which you have the executable and run the command
```bash
gprof -l nonlinearsigma gmon.out > profiling_results.txt
```
You will then be able to see the profiling report. In general, you should set "USE_GPROF" to false, unless you are looking to optimize the code or troubleshoot it.

The flag "USE_TEST_PRINT_STATEMENTS" activates print statements throughout the code. This is useful for debugging, but should generally be set to FALSE as it slows down the code.

If you run into major problems, set "USE_EXTREME_TEST_CONDITION" to TRUE. This will run a testing suite built into the code, but will not run the usual simulation. This can help you identify problems in the code, and the testing suite is a function inside the main function, which can be modified as needed to add more tests.

To switch from using arcsin to calculate $Q_{L}$ to using arccos, set "USE_CHECK_QL_COS" to TRUE. In general, this should be set to FALSE< as we want to use arcsin due to its useful symmetry.

Finally, if you want to remove the random number generation and use a constant value for the random numbers, set "USE_CONST_RN" to true. This should only be done when testing the code.

#### Compiling the code

Once you have set the flags you want to use, you can compile the code by typing

```bash
make -f make_sigma
```
It can be useful to run 

```bash
make -f make_sigma clean
```
first to clear out old .o files that might not be updated otherwise.

Once you have compiled the code, you should have an executable by the name of ```nonlinearsigma``` which you will use to run the simulation.

#### Other important things to note before running the code

The simulation requires a list of parameters. We give those to the code through the text file ```inputs.txt```. That file consists of a list of parameter keywords and a value. Be very careful not to change the format of this document, only the numbers, otherwise the simulation won't be able to understand what values get assigned to what parameters. 

The file ```inputs.txt``` should look like this

```
L = 10
beta = 1.6
itheta = 1.
ntherm = 1000
nMC = 1000
freq = 100
```
Where ```L``` gives the length of the square lattice, ```beta``` is $\beta = 1/g_{L}$ and should be set to 1.6, ```itheta``` is the imaginary value given for the topological term and is given in fractions of $\pi$ (so you should enter 0.5 if you want $i \theta = \pi/2$), ```ntherm``` is the number of steps you want the simulation to take for thermalization, ```nMC``` is the number of steps in the Monte Carlo loop after thermalization, and ```freq``` sets the number of steps between saved configurations.

#### Running the code -- single job in interactive node

If you want to run this in an interactive node to test, you can request an interactive node with the following command:

```bash
salloc --cpus-per-task=1 --time=00:30:00
```
where you can adjust the number of cpus per task if you want to run the parallelized code and the time requested is in format hh:mm:ss

You can run the job using the command

```bash
./nonlinearsigma inputs.txt
```
```nonlinearsigma``` is the name of the executable created when you compiled the code, and inputs.txt is the list of inputs mentioned above. 

#### Running the code -- single SLURM submission 

Once you're comfortable with the code and want to submit a job that's longer than your interactive session (or one that may need to run overnight, etc), you do that by writing a SLURM script and sending that script to the scheduler.

The script is called ```submit_sigma.sh``` and looks like this:
```bash
#SBATCH --job-name=nonlinearsigma_omp_test           # Job name
#SBATCH --mail-type=ALL                              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cberger@smith.edu                # Where to send mail
#SBATCH --partition=phyq                             # Which partition to use
#SBATCH --nodes=1                                    # Number of nodes
#SBATCH --cpus-per-task=1                           # Number of threads per task (OpenMP)
#SBATCH --mem=1gb                                    # Job memory request
##SBATCH --time=05:00:00                             # Time limit hrs:min:sec
#SBATCH --output=nonlinearsigma_omp_test_%j.log      # Standard output 
#SBATCH --error=err_nonlinearsigma_omp_test_%j.log   # Standard output and error log

pwd; hostname; date

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Running nonlinear sigma on single CPU core"

/usr/bin/time -v ./nonlinearsigma inputs.txt

date
```

Let's walk through what this script does.

The first line:
```bash
#SBATCH --job-name=nonlinearsigma_omp_test           # Job name
```
Assigns a name to the job, which can help you keep track of what is running in the queue. I tend to name all my jobs in a certain phase of the work the same thing (e.g. "nonlinearsigma_small_L_tests" if I'm testing out the script on small lattices or "nonlinearsigma_first_production_run" if I'm starting to take data for real). This is just for your own information, so name it whatever you'd like.

The next two lines
```bash
#SBATCH --mail-type=ALL                              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cberger@smith.edu                # Where to send mail
```
tell SLURM who to email and when to email you. I have it set to email me at my Smith email anytime a job starts, ends, or fails. I find this helpful, especially if trouble arises, but it's your choice what to put here (just don't leave my email address in!).

The next line
```bash
#SBATCH --partition=phyq                             # Which partition to use
```
sends the jobs to the partition that belongs to the physics department. We have priority on this node, but we can request time on other nodes if we really need to. We'd need to talk to CATS about that if we wanted to do it.

The next lines specify the code's needs:
```bash
#SBATCH --nodes=1                                    # Number of nodes
#SBATCH --cpus-per-task=1                           # Number of threads per task (OpenMP)
#SBATCH --mem=1gb                                    # Job memory request
```
Since we are using OpenMP for parallelization, we only need one node, but we will want to change ```bash --cpus-per-task``` to something larger for OpenMP. Unless doing a very large lattice (L > 100), I tend to ask for 30 CPUs per task, which allows for 2 jobs to run at a time on each node. If you're doing something large, you might want to consider asking for 60 CPUs per task, which will occupy an entire node.

This line
```bash
##SBATCH --time=05:00:00                             # Time limit hrs:min:sec
```
sets a time limit -- it will cut off your code when that limit is reached, whether it is done or not. Note there is an extra `#` here -- that means it's commented out, so it will not have a time limit. I tend to only use the time limits when testing.

These next two lines
```bash
#SBATCH --output=nonlinearsigma_omp_test_%j.log      # Standard output 
#SBATCH --error=err_nonlinearsigma_omp_test_%j.log   # Standard output and error log
``` 
give names to the output and error files. These files are created when the code runs and can help you debug if things go wrong.

All the above were instructions for SLURM, which schedules jobs on the machine's available resources. Now we get into the commands to actually run the code.

```bash
pwd; hostname; date
```
This just prints where the job is being run from and the date and time before the job starts.

```bash
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
```
this sets the number of threads in OpenMP equal to the number of CPUs per task you assigned above.

```bash
echo "Running nonlinear sigma"

/usr/bin/time -v ./nonlinearsigma inputs.txt
```
This prints out an annoucement that you're starting the job, which is useful in the logfile, and then runs the simulation with the inputs file. Make sure you have the inputs file and the executable in the same folder as the SLURM script so the computer can find them.

And finally
```bash
date
```
we print the datetime stamp at the end of the simulation.

#### Running the code -- batch SLURM submissions

Each inputs file is one set of parameters, and we need to get lots of data. If you're not interested in manually setting up these SLURM scripts, input files, etc by hand, I don't blame you. That's why I wrote a Python code to do most of the work for you.

In this code, you specify what parameters you want to run, and the script creates all the appropriate directories and puts the correct inputs file, slurm script, and the executable in each directory. You still have to go in and submit the files yourself, but it's much easier than writing all these scripts yourself.

The script is called ```create_input_files.py``` and here is the part you will need to modify:

```python
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
itheta_list = [0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.,1.125]

script_name = "nonlinearsigma"
job_name = "nlsigma_prelim_tests"
email = "cberger@smith.edu"
num_cpus = 30

```
Beta, ntherm, nMC, and freq all take one number as input, but you can create a list of the number of lattice lengths you want, and the values for itheta (remember these are fractions of pi). 

```script_name``` is the executable, so it should be "nonlinearsigma", but ```job_name``` is what will be put in for the job name in the SLURM script. Similarly, you can enter the email address you want included in the SLURM script and choose how many CPUs you want.

Once you've modified this script to have the values you want, put it in a directory for this batch. I tend to name those something like ```run_yyymmdd```. Also in that directory should be the executable, so copy that in once you've compiled the code.

Then, inside the batch directory, go ahead and run the python script

```bash
python create_input_files.py
```

When it's done, you should see the subdirectories created -- one for each job. You need to go into each subdirectory to submit the jobs using

```bash
sbatch submit_sigma.sh
```
and it will submit the jobs to the scheduler.

You can run this with 
```bash
sbatch submit_sigma.sh
```
and you can check the status of your jobs any time with the command
```bash
squeue -u your_username
```

## Data and Analysis


```python

```

## Current issues and open questions


```python

```
