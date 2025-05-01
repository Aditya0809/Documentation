# Software Installation

## Petsc

* Clone the [Petsc](https://petsc.org/release/install/download/) repository
* Obtain the modified files from the "changes-petsc" directory in this repository
* In the cloned PETSc directory, replace euler.c with the modified version from "changes-petsc". This file is located at /path-to-petsc/src/ts/impls/explicit/.
* Replace tsimpl.h similarly, with its modified version found in "changes-petsc", located at /path-to-petsc-3.15.4/include/petsc/private/tsimpl.h/.
* Follow the configuration and compilation instructions provided at: [petsc_installation](https://petsc.org/release/install/install/.) Configuration options may vary; for example, I use: OPTFLAGS="-Ofast -march=native -mavx2 -mfma -fno-finite-math-only" ./configure --with-debugging=no COPTFLAGS="$OPTFLAGS" CXXOPTFLAGS="$OPTFLAGS" FOPTFLAGS="$OPTFLAGS" --with-blaslapack-dir=$MKLROOT --with-scalapack-lib="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64" --with-scalapack-include=$MKLROOT/include
* Set the PETSC_DIR and PETSC_ARCH environment variables to their appropriate paths (visible upon successful configuration).

## PetIGA

This project also leverages a specialized version of PetIGA, initially developed by L. Dalcin et al. . To install it, follow these directions:

* Navigate to the mod_PetIGA folder you've downloaded.
* Execute "make all" and then "make test" for compilation and testing.
* Define the PETIGA_DIR and PETIGA_ARCH environment variables with their corresponding paths.