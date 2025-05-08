# Software Installation

Here is a detailed guide for installing PETSc, modified PetIGA, and IgaKit on the Bridges2 cluster.

## Necessary modules

1.  mkl
2.  intelmpi/2021.3.0-intel2021.3.0
3.  anaconda (for IGAKIT) 


## PETSc

1. **Clone the PETSc repository**  : [Download PETSc](https://petsc.org/release/install/download/)

2. **Copy modified files from this repository**  
   This project uses a modified version of PETSc to support **explicit time stepping with a lumped mass matrix** in the IGA context.

    To enable this, two core PETSc source files have been updated:

    - `euler.c` : located in `src/ts/impls/explicit/`
    - `tsimpl.h`: located in `include/petsc/private/`

    You can download the modified versions from this repository: [changes_petsc](files/changes_petsc.zip)

    🔧 **Summary of Changes**

    ✅ `tsimpl.h` (structure update)

    A new vector field for the lumped mass matrix has been added to the internal `TS` structure:

    ```c
    /*----------------------- Lumped Mass Matrix -----------------------------------*/
    Vec vec_lump;
    ```
    This line is inserted inside the `_p_TS` structure to store the vector of lumped mass coefficients. 

    ✅ `euler.c` (RHS update in TSStep_Euler())

    The right-hand side (RHS) vector is now divided by the lumped mass vector (`vec_lump`) before updating the solution, as shown below:

    ```c hl_lines="2 3 10 11 12"
        TS_Euler       *euler = (TS_Euler*)ts->data;
        Vec            solution = ts->vec_sol, update = euler->update;
        Vec            lmc = ts->vec_lump, update2 = euler->update;
        PetscBool      stageok, accept = PETSC_TRUE;
        PetscReal      next_time_step = ts->time_step;
        PetscErrorCode ierr;

        PetscFunctionBegin;
        ierr = TSPreStage(ts, ts->ptime);CHKERRQ(ierr);
        ierr = TSComputeRHSFunction(ts, ts->ptime, solution, update);CHKERRQ(ierr);
        ierr = VecPointwiseDivide(update2, update, lmc);
        ierr = VecAXPY(solution, ts->time_step, update2);CHKERRQ(ierr);
        ierr = TSPostStage(ts, ts->ptime, 0, &solution);CHKERRQ(ierr);
        ierr = TSAdaptCheckStage(ts->adapt, ts, ts->ptime, solution, &stageok);CHKERRQ(ierr);
        if (!stageok) {
        ts->reason = TS_DIVERGED_STEP_REJECTED;
        PetscFunctionReturn(0);
        }
        ierr = TSFunctionDomainError(ts, ts->ptime + ts->time_step, update, &stageok);CHKERRQ(ierr);
        if (!stageok) {
        ts->reason = TS_DIVERGED_STEP_REJECTED;
        PetscFunctionReturn(0);
        }
        ierr = TSAdaptChoose(ts->adapt, ts, ts->time_step, NULL, &next_time_step, &accept);CHKERRQ(ierr);
        if (!accept) {
        ts->reason = TS_DIVERGED_STEP_REJECTED;
        PetscFunctionReturn(0);
        }
        ierr = VecCopy(solution, update);CHKERRQ(ierr);
        ierr = VecCopy(solution, ts->vec_sol);CHKERRQ(ierr);

        ts->ptime += ts->time_step;
        ts->time_step = next_time_step;
        PetscFunctionReturn(0);
    ```

    These changes enable the use of lumped mass matrix-based explicit integration schemes, which are particularly useful when using PetIGA for transient simulations.




3. **Replace PETSc source files**  
   - Change `tsimpl.h` file located at:
     ```
     /path-to-petsc/include/petsc/private/tsimpl.h
     ```
   - Overwite the original euler.c found at:
     ```
     /path-to-petsc/src/ts/impls/explicit/euler.c
     ```


4. **Configure and compile PETSc**  
   Follow the official guide: [PETSc Installation](https://petsc.org/release/install/install/)  
   Example configuration:
   ```bash
   OPTFLAGS="-Ofast -march=native -mavx2 -mfma -fno-finite-math-only"
   ./configure --with-debugging=no \
       COPTFLAGS="$OPTFLAGS" \
       CXXOPTFLAGS="$OPTFLAGS" \
       FOPTFLAGS="$OPTFLAGS" \
       --with-blaslapack-dir=$MKLROOT \
       --with-scalapack-lib="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64" \
       --with-scalapack-include=$MKLROOT/include
   make all
   ```

5.  **Set environment variables**
    After successful configuration, define the following:
    ```bash
        export PETSC_DIR=/path/to/petsc
        export PETSC_ARCH=your-architecture-name
    ```




## PetIGA

This project uses a modified version of **PetIGA**, originally developed by L. Dalcin et al. You can download our version here: [Download mod_PetIGA.zip](files/mod_PetIGA.zip) 

1.  Navigate to the `mod_PetIGA` folder you've downloaded. 
2.  Run: ```make all make test ```  for compilation and testing
3.  Set the following environment variables: 
    ```bash
        export PETIGA_DIR=/path/to/mod_PetIGA 
        export PETIGA_ARCH=your-arch 
    ``` 


### Why we modified PetIGA 

The original PetIGA implementation performs explicit time stepping as: \\[ U^{n+1} = U^n + \Delta t \cdot \mathcal{R}(U^n) \\] 
However, this ignores the presence of the **mass matrix**, which should appear in the weak form of time discretization. 
The proper discretized formulation is: \\[ M \cdot U^{n+1} = M \cdot U^n + \Delta t \cdot \mathcal{R}(U^n) \\] 
To avoid the cost of inverting \\( M \\), we apply the **lumped mass matrix** technique, where: 

\[
\mathcal{M}_{AB} = 0  \text{if} A \neq B 
\]


This converts the system into a diagonal form, allowing for efficient inversion: \\[ U^{n+1} = U^n + \Delta t \cdot \mathcal{M}^{-1} \mathcal{R}(U^n) \\] 
This modification improves performance in explicit schemes while maintaining physical correctness.

\[
    \mathcal{M}_{AB}=\begin{cases}
          \sum_{b=1}^{n_b} {M}_{Ab} \: & \text{if} \; A=B, \\
          0 \: & \text{if} \; A\neq B. \ 
     \end{cases}
\]


\[ 
 \mathcal{M}_{AB} = \sum_{b=1}^{n_b} {M}_{Ab} \,\, \text{if} \,\, A=B,
\]

### Extra note on the **explicit RHS helper**

In the stock PetIGA API you call `IGACreateTS()` (implicit, Jacobian‑based) or `IGACreateTS2()` (matrix‑free).  
For purely explicit schemes we added a convenience wrapper:

```c
    /* Purpose: create a PETSc TS object, attach the IGA,
    *          and wire the **right‑hand side only** callback
    *          required for explicit time steppers (TSEULER, TSRK, etc.) */
    PetscErrorCode IGACreateTS3(IGA iga, TS *ts)
    {
    MPI_Comm       comm;
    Vec            U;
    Vec            F;
    Mat            J;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
    PetscValidPointer(ts,2);

    ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
    ierr = TSCreate(comm,ts);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)*ts,"IGA",(PetscObject)iga);CHKERRQ(ierr);
    ierr = IGASetOptionsHandlerTS(*ts);CHKERRQ(ierr);

    ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
    ierr = TSSetSolution(*ts,U);CHKERRQ(ierr);
    ierr = VecDestroy(&U);CHKERRQ(ierr);

    ierr = IGACreateVec(iga,&F);CHKERRQ(ierr);
    ierr = TSSetRHSFunction(*ts,F,IGATSFormRHSFunction,iga);CHKERRQ(ierr);
    ierr = VecDestroy(&F);CHKERRQ(ierr);

    PetscFunctionReturn(0);
    }
```

Key steps inside `IGACreateTS3()` :

| Step | What it does |
|------|--------------|
| **1** | `TSCreate` and attach the IGA object via `PetscObjectCompose`. |
| **2** | Allocate a global vector **U** and register it with `TSSetSolution`. |
| **3** | Allocate a work vector **F** and register **only the RHS** via `TSSetRHSFunction`; no Jacobian is set. |

The actual residual assembly is delegated to `IGATSFormRHSFunction()`, so later you simply call:

```c
    TSSetType(ts, TSEULER);   /* or TSRK, etc. */
    TSSolve(ts, U);           /* PETSc advances with explicit lumped‑mass update */
```



## IGAKit for Visualization

To download and install IGAKit, follow the instructions provided in the [official GitHub repository](https://github.com/dalcinl/igakit).

**IGAKit** is a lightweight Python toolkit that complements PetIGA by providing:

* **Geometry generation & manipulation**  
  – construct NURBS curves, surfaces, and volumes directly in Python.

* **Visualisation helpers**  
  – export control meshes, knot lines, and solution fields (VTK) for ParaView or VisIt.

* **I/O utilities**  
  – read/write `.iga` and IGES files, perform uniform refinements, and inspect knot vectors.

In this documentation we will use IGAKit mainly to **visualise simulation output** produced by PetIGA by converting `.dat` files to `.vtk` files, which can be visualized in paraview.  


**Installation tip:** IGAKit is pure‑Python. Activate your conda environment on Bridges2 and run  
```bash
     pip install --user igakit
```  
or clone the repo and install manually:  
```bash
     python setup.py install --user
```