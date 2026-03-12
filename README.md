# 2D Lattice Boltzmann Method Flow Simulation (MATLAB)

This repository provides a compact MATLAB implementation of a **2D Lattice Boltzmann Method (LBM)** solver using the **D2Q9 lattice model** and **BGK collision operator**.

The code simulates **pressure-driven channel flow past a circular obstruction**, producing the classical **Kármán vortex street**.

The goal of this project is to help people come into LBM world with a simple MATLAB script: with only a few hundred lines of MATLAB code, users can run a complete fluid simulation and visualize velocity, pressure, and vorticity fields in real time.

---

# How to Run

1. Download the .m code.

2. Open MATLAB

3. Run the .m script


The simulation will automatically:

- initialize the computational domain
- run the LBM solver
- visualize velocity, pressure, and vorticity fields
- optionally export GIF animations of the flow evolution

---

# Full Documentation

A detailed explanation of the **algorithm, equations, and code structure** is available here:

👉**Project Documentation**

https://zongze-li-pon.github.io/completed-projects/2019-11-14-2d-lattice-boltzmann-method-flow-simulation-matlab/

The document includes the explanation of all code implementation:

- LBM governing equations
- D2Q9 lattice model
- pressure boundary implementation
- explanation of each code section
- physical interpretation of the simulation
