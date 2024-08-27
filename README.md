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
The -DDATA flag will tell the program to save all the required data, which, at the moment, consists only of the MCS, the magnetization and the Binder cummulant.

```
#ifdef DATA
	//fprintf(arq1, "%d,%f,%f,%f,%f,%f,%f,%f,%f\n", mcs, ET, M, M2, M4, U, ET2, CV, EE);
	fprintf(arq1, "%d,%f,%f\n", mcs, M, U);
#endif
```
