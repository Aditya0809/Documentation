# Welcome to the Explicit PetIGA Documentation

[PetIGA](https://github.com/dalcinl/PetIGA) extends the [PETSc](https://petsc.org/release/) toolkit by adding NURBS‐based isogeometric discretisation and automated form integration.  With it you can build 2‑D and 3‑D solvers—linear or nonlinear, steady or transient—without writing domain‑decomposition, MPI halo‑exchange, or low‑level assembly code.  The major aim while using PetIGA is to code the residual (or right‑hand side in explicit cases) and, if needed, the Jacobian at each Gauss point; PetIGA and PETSc handles everything else.



## Advantages of explicit time stepping 

| Benefit | Comment |
|---------|---------|
| **Low memory footprint** | No global matrix inversion; ideal when the number of degrees of freedom (DOFs) is very large. |
| **No Jacobian required** | Saves coding effort and sidesteps numerical differentiation costs. |
| **Embarrassingly parallel** | Each time step is mostly local, giving excellent scalability. |
| **Simple to implement** | Update rule is a vector operation; easier for newcomers. |

> **Caveat:** stability is dictated by a CFL‑type time‐step limit.



## PetIGA’s capabilities

* **NURBS / B‑spline basis** for exact CAD geometry and high‑order, \(C^{p-1}\) continuity.
* **Unified API** for linear, nonlinear, and time‑dependent PDEs.
* **Strong scaling demonstrated** to 2048 CPU cores on problems in solid & fluid mechanics.
* **Automatic numerical differentiation** available when an analytic Jacobian is tedious to formulate.
* **Open‑source, MIT licence** – easy to comprehend, extend, and reuse alongside PETSc’s solvers and preconditioners.



## What you’ll find in this documentation

1. **Installation guide** - install and compile the modified PetIGA and Petsc
1. **Quick‑start tutorials** – running your first 2‑D heat conduction example.  
3. **Worked examples** in fourth-order Cahn Hilliard equation.  
4. **Scaling tips** – how to push your model to hundreds of millions of DOFs.  
5. **Reference pages** – installation flags, command‑line options, and helper scripts.

Whether you are a graduate student prototyping a new constitutive law or an engineer scaling to thousands of cores, this guide will get you productive with PetIGA and PETSc.

