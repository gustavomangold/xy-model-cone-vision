# XY model with cone vision

This repository consists of code from the paper _Long-range Order and Directional Defect Propagation in the Nonreciprocal XY
Model with Vision Cone Interactions_, Loos et. al..

## Python implementation

We first implement the model in python in a more explicit manner, obtaining results from both the pure metropolis algorithm as well as using Glauber dynamics, which is the method from the paper.
The python program outputs only magnetizations, as it is more for testing the implementation.

## C implementation

This version is a lot faster, and should be used for running the actual simulations. 
