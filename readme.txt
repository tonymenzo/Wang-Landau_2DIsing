Final Project - Wang Landau Sampling of the 2D Ising Model

Programmer: Tony Menzo
Last Modified: 4/27/2018

Summary: This program takes a spin 1/2 two-dimensional lattice and implements the Wang-Landau(WL) sampling algorithm to determine the density of states "function" for the system. This algorithm is an alternative to the popular Metropolis algorithm. The main difference between the two is that in WL sampling the random walk is performed in energy space. 
Once this function is determined it is reletively straightforward to calculate the partition function of the ensemble along with other thermodynamic variables such as internal energy, Cv (specific heat), entropy, and Gibbs free energy.


Welcome!!!! You are about to explore my program implementing WL sampling!

Here is a guide

(1) Open up "ising_model_WL.cpp" 
(2) Compile and link the program using "make -f make_ising_model_WL"
(3) I've set the side length L to 15 to start, this tells the program how big of a lattice you want, start out from 15 and go up and down from there, around L=16 the program begins to have larger run times.
(4) My main focus for this project was to implement the algorithm and obtain the density of states function, I was successful in obtaining g(E). 
(5) It was suggested that the intial f be set to 2.71 and terminate at 1, this gives enough interations in order to have enough data for g(E)
 

How the program works: 
(1) The program creates a configuration with all of the spins up to calculate the maximum energy of lattice with side length L 
(2)Then a random configuration is generated and the energy of that configuration is set as the current configuration. This energy and configuration will be be the first energy used in our monte carlo simulation.
(3) Initialize both counting arrays, g[E] and hist[E]. We are using ln(g[E]) and ln(f) in order to avoid double precision round off error that arises from many multiplications. See source (2). The modifcation factor is is reduced by a power of 1/2 after each iteration of the monte carlo simulation
(4) Now we start the Monte Carlo Simulation, we do 10,000 simulations and then check for the flatness of the histogram. We do this until the modification factor is <= exp(10^-8)
(5) After the Monte Carlo simulation we need to normalize g(E)
(6) Finally, I interpolate the the data.



Files: 
"wanglandau.c" - this code written in C implementing the Wang-Landau algorithm that I used partially as a reference. My main resource was "A Survey of Computational Physics"

"random_seed.cpp" - this generates a random seed number that the gsl uniform random number generator uses to generate it's random number. 

"partition.dat" - this is where data for the partition function is input (still working on this part of the code)

"make_ising_model_WL" - this is the makefile for the program

"L_time.dat" - this file shows how the runtime scales with size of the lattice

"ising_model_WL.cpp" - this is where the magic happens

"GslSpline.cpp/GslSpline.h" - this files contain the the class and class headers for the splining routines I used when fitting g(E) (and eventually Z(T))


"g_E_plt.ps" - this shows my plots for the density of states for various L

"g_E_histogram_plt.ps" - this shows the total histogram illustrating the number of times each energy level was visited. Flat histogram means my program worked!! All of the energy levels were visited equally. 

"g_E_histogram.plt/g_E.plt" These are my plot files you can load in gnuplot to see g(E) and the histogram

"g_E_fit.dat" - this is my data from the interpolation of my calculated g(E)

"E_vs_kT.plt/.ps" - these files are to illustrate how energy scales with temperature (still working on this part) 

Upgrades at a later date:

There are a few things I am still going to work on for this routine. I plan on calculating all of the thermodynamic variables, as well as fitting the data rather than interpolating. 
Some of the routines can be optimized much like those in ising_opt.cpp. 
I also think in it's current state is pretty long, much of the routines like the monte carlo simulation could be implemented in a class. 
Things I'd like ot put into a class:
Monte Carlo Simulation
Printing data
Fitting 
I would also like to generalize the program to include 1D and 3D lattices. 

Thank you for for your help and guidance this semester. This has been my favorite course of my undergraduate career!!!

 
