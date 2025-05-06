# Homepage

[PetIGA](https://github.com/dalcinl/PetIGA) is a software framework that approximates the solution of partial differential equations using isogeometric analysis. It is an extension of PETSc, adding the NURBS discretization capability and the integration of forms. This framework can be used to solve linear, nonlinear, time-dependent, or time-dependent nonlinear problems. In this framework, the user has to provide the evaluation of the linear form (right-hand side, or residual of a nonlinear problem) at a Gauss point, as well as the bilinear form (left-hand side or Jacobian of the nonlinear residual). The framework is designed so that researchers can focus on the physics of the problem and ignore issues of parallelism and performance.


# Welcome to the PetIGA + PETSc Documentation

**PetIGA** extends the PETSc toolkit by adding NURBS‐based isogeometric discretisation and automated form integration.  With it you can build 2‑D and 3‑D solvers—linear or nonlinear, steady or transient—without writing domain‑decomposition, MPI halo‑exchange, or low‑level assembly code.  Your task is reduced to coding the residual (or right‑hand side) and, if needed, the Jacobian at each Gauss point; PetIGA and PETSc handle everything else.

---

## Why this guide focuses on **explicit** schemes first  

| Benefit | Comment |
|---------|---------|
| **Low memory footprint** | No global matrix inversion; ideal when the number of DOFs is very large. |
| **No Jacobian required** | Saves coding effort and sidesteps numerical differentiation costs. |
| **Embarrassingly parallel** | Each time step is mostly local, giving excellent scalability. |
| **Simple to implement** | Update rule is a vector operation; easier for newcomers. |

> **Caveat:** stability is dictated by a CFL‑type time‐step limit.  We’ll show you how to choose a safe Δt.

---

## Snapshot of PetIGA’s capabilities

* **NURBS / B‑spline basis** for exact CAD geometry and high‑order, \(C^{p-1}\) continuity.
* **Unified API** for linear, nonlinear, and time‑dependent PDEs.
* **Strong scaling demonstrated** to 4 096 CPU cores on problems in solid & fluid mechanics.
* **Automatic numerical differentiation** available when an analytic Jacobian is tedious.
* **Open‑source, MIT licence** – easy to hack, extend, and reuse alongside PETSc’s solvers and preconditioners.

---

## What you’ll find in this documentation

1. **Quick‑start tutorials** – from compiling PETSc + PetIGA to running your first 2‑D diffusion example.  
2. **Explicit time‑stepping recipes** – choosing Δt, managing output, and monitoring stability.  
3. **Worked examples** in solid (elastic plate) and fluid (2‑D advection–diffusion) mechanics.  
4. **Scaling tips** – how to push your model to hundreds of millions of DOFs.  
5. **Reference pages** – installation flags, command‑line options, and helper scripts.

Whether you are a graduate student prototyping a new constitutive law or an engineer scaling to thousands of cores, this guide will get you productive with PetIGA and PETSc.

#### Code for a specific language

Some more code with the `py` at the start:

``` py
import tensorflow as tf
def whatever()
```

#### With a title

``` py title="bubble_sort.py"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

#### With line numbers

``` py linenums="1"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

#### Highlighting lines

``` py hl_lines="2 3"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

## Icons and Emojs

:smile: 

:fontawesome-regular-face-laugh-wink:

:fontawesome-brands-twitter:{ .twitter }

:octicons-heart-fill-24:{ .heart }