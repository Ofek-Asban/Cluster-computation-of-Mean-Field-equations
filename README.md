# Cluster-computation-of-Mean-Field-equations
## Physics Background:
The project deals with the time evolution (relaxation) of the Electron-glass (EG) and dipolar Spin-glass (SG) systems. The EG depicts a system of electrons situated in random locations with a random on site energy. The electrons interacts via the couplomb interaction. The dynamics ofthe system come from quantu, tunneling of the electrons between sites asysted by phonons (matter waves). 
The dipolar Spin-glass system has similar strucure only instead electrons we have spins in each site and theor dynamics comes from tunneling between down and up orientations within the same site. For further details see the paper [Effect Of Interactions and Disorder](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.144207)

To find the thermal (or ground state) the mean field approximation is applies. In this approximation the problem reduces to single particle problen without interactions. The some comtribution of the interations is embedded in the effective on-site energy and can be found by the solution of highly non-linear self consist set of equations known as the Mean-Field equations. 

## Goal of the Algorithm:
The goal is to find the solution of the self-consistent equations of the EG and SG systems.

The EG Mean-Field equations are, 
$$(1) \\ E_i = \epsilon_i + \frac{1}{2}\sum_{j (\neq j)} \frac{e^2}{r_{ij}}(e^{E_j/T} + 1)^{-1}.$$
Similarly the SG Mean-Field equations are,
$$(2) \\ \Delta_i' = \Delta_i + \frac{1}{4}\sum_{j (\neq j)} \frac{u_{ij}}{r^3_{ij}}\tanh\left(\frac{\Delta_j'}{2T}\right),$$
where the on-site energy $\epsilon_i$ and the spin bias energy $\Delta_i$ are modeled as uniformly distributed numbers, $r_{ij}$ is the distance between sites i and j in two dimentions and thus represented by symmetric and uniformly distributed matrix with 0 diagonal (with periodic boundary condisitons), and $u_{ij}$ is a Gaussian distributed parameter with controlled variance.

## Structure of the Algorithm
One of the methods to solve equations (1) and (2) (i.e. find values of $E_i$ and $\Delta_i'$ that satisfy equations (1) and (2)). I chose to use relaxation method which contains two steps:

(a) Initializing $E_i$ and $\Delta_i'$ with uniform random values and iterating back and forth between the left and right hand side of equations (1) and (2) untill the difference between iterations is small enough.

(b) To increase the the accuracy of the solution once a solution is found the systems is heated and and cooled again (increasing and then lowering the value of T variable to its desired value). This allows the system to overcome large barriers and find a better minimum point (the energy landscape has many local minima).

Since the systems has random variables we our algorithm consists of a third step 
(c) Perform steps (a) and (b) many realizations of the random variables in the model and then average over the solutions $E_i$ and $\Delta_i'$. The final representation is a normalized histogram (distribution) of $E_i$ and $\Delta_i'$ vectors presented in 

The three steps are computationally costly, comutation time of steps (a) and (b) alone can take between roughly 5-30 minutes on a single desktop computer (depending on the random values of the relization, and the initialization values of $E_i$ and $\Delta_i'$). Thus, to reduce the computation time from days to hours I used distributed computing in SGE cluster and sent 100 to 100,000 of random realizations to be processed in parallel and then made the average, using scripts that automated the process. 

The number of realizations, n, is chosen such that $n \times N = 10^6$ where N is the the number of sites in the EG system if spins in the SG system. 

## Files in the repository

The files presented here are for the solution of Eq (2). The same procedure is done for solving Eq (1).

