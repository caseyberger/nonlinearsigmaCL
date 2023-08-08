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

Once you have compiled the code, you should have an executable by the name of ```nonlinear_sigma``` which you will use to run the simulation.

#### Running the code -- single job in interactive node

#### Running the code -- SLURM submission and batching


```python

```
