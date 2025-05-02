# Software Installation

## PETSc

1. **Clone the PETSc repository**  
   [Download PETSc](https://petsc.org/release/install/download/)

2. **Copy modified files from this repository**  
   This project uses a modified version of PETSc to support **explicit time stepping with a lumped mass matrix** in the IGA context.

    To enable this, two core PETSc source files have been updated:

    - `euler.c` : located in `src/ts/impls/explicit/`
    - `tsimpl.h`: located in `include/petsc/private/`

    You can download the modified versions from this repository: [changes_petsc](files/changes_petsc.zip)

### 🔧 Summary of Changes

#### ✅ `tsimpl.h` (structure update)

    A new vector field for the lumped mass matrix has been added to the internal `TS` structure:

    ```c
    /*----------------------- Lumped Mass Matrix -----------------------------------*/
    Vec vec_lump;
    ```
   - This line is inserted inside the `_p_TS` structure to store the vector of lumped mass coefficients. 

#### ✅ `euler.c` (RHS update in TSStep_Euler())

 The right-hand side (RHS) vector is now divided by the lumped mass vector (`vec_lump`) before updating the solution, as shown below:

  ```c
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
   - Download the modified `euler.c` file and overwrite the original:
     ```
     /path-to-petsc/src/ts/impls/explicit/euler.c
     ```
   - Download the modified `tsimpl.h` file:
     ```
     /path-to-petsc/include/petsc/private/tsimpl.h
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

This project uses a modified version of PetIGA, initially developed by L. Dalcin et al. Follow the steps below to install it: [Download](files/mod_PetIGA.zip)

1. Navigate to the mod_PetIGA folder you've downloaded.
2. Execute "make all" and then "make test" for compilation and testing.
3. Define the PETIGA_DIR and PETIGA_ARCH environment variables with their corresponding paths.

