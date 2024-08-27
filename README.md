# XY model with cone vision - Mini Documentation

This repository consists of code from the paper _Long-range Order and Directional Defect Propagation in the Nonreciprocal XY
Model with Vision Cone Interactions_, Loos et. al..

## Python implementation

We first implement the model in python in a more explicit manner, obtaining results from both the pure metropolis algorithm as well as using Glauber dynamics, which is the method from the paper.
The python program outputs only magnetizations, as it is more for testing the implementation.

## C implementation

This version is a lot faster, and should be used for running the actual simulations. 
The implementation consists of two different files, one legacy header _mc.h_ that has Monte Carlo functions and the _xycone.c_ file, which implements the actual simulation.

Compilation should be done with gcc, using the command _gcc -DDATA -Wall c-version/xy_cone.c -lm -O3_ and then ran with ./a.out TEMP THETA (degrees) SEED.
The -DDATA flag will tell the program to save all the required data for the _TMAX_ steps (which will run after the transient steps, defined as _TRAN_). In the code, it currently consists only of the MCS, the magnetization and the Binder cummulant, but can save anything that's updated throughout runtime.

```
#ifdef DATA
	//fprintf(arq1, "%d,%f,%f,%f,%f,%f,%f,%f,%f\n", mcs, ET, M, M2, M4, U, ET2, CV, EE);
	fprintf(arq1, "%d,%f,%f\n", mcs, M, U);
#endif
```
The process of simulation occurs as it is defined in the ```main()``` function. 
After allocating space for each spin and initializing them randomly, as well defining its neighbours inside the function  ``` void initialize ```, we then loop through all transiente Monte Carlo steps, iterating inside ```void mc_routine```, which in turn samples _L*L_ spins.
Each spin sampled has its energy _Ei_ calculated, and then a new value is randomly assigned to its orientation, with:
```  
flip = (spin[i] + (2*M_PI)*FRANDOM);

int_flip = (360 * (flip / (2*M_PI)));

int_flip = int_flip % 360;

flip     = (int_flip * 1.0 * (2*M_PI/360));
```

In the first line the random spin is sampled uniformly, and the three extra steps are done in order to rotate it back to the (0, 2pi] interval, in order to maintain consistency and prevent future work when printing states and comparing angles for the cone vision inequalities.

Then, we calculate the energy _Ef_ of the new spin and accept it with probability _G_, given by:
```
dE = Ei-Ef;

G = (1-tanh(dE/(2*TEMP)))/2;

if(FRANDOM < G)
{
	spin[i] = flip;
}
```
After all _L*L_ samples, 
