# Software Installation

## PETSc

1. **Clone the PETSc repository**  
   [Download PETSc](https://petsc.org/release/install/download/)

2. **Copy modified files from this repository**  
   Get the files from the `changes-petsc/` directory in this repository.

3. **Replace PETSc source files**  
   - Download the modified `euler.c` file and overwrite the original: [Download euler.c](files/euler.c)  
     ```
     /path-to-petsc/src/ts/impls/explicit/euler.c
     ```
   - Download the modified `tsimpl.h` file: [Download tsimpl.h](files/tsimpl.h), and overwite [Download](files/changes-petsc.zip)
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

