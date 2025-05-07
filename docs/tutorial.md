# Tutorial I · 2‑D Heat‑Conduction Example

**Objective** – show an *explicit* PetIGA implementation for the transient heat‑diffusion equation and how to run it on your cluster.

---
## 1 · Problem description

| Item | Value / comment |
|------|-----------------|
| **Domain** | 2‑D square, side length \(L = 50\). |
| **Boundary conditions** | *Periodic* on all four edges (model a repeating tile). |
| **Initial temperature** | Parabolic “hot spot” centred at \((L/2,L/2)\):<br> \(T(\mathbf x,0)=100\Bigl[\,1-\bigl(\tfrac{r}{0.5L\sqrt{2}}\bigr)^2\Bigr]\) where \(r=\|\mathbf x-\mathbf x_c\|\) and values below 0 are clipped. |
| **Material property** | Thermal diffusivity \( \alpha = \dfrac{\kappa}{\rho c_p}=5.0 \). |
| **Discretisation** |  \(\Delta x = 0.5\) ⇒ \(100\times100\) elements; time step \(\Delta t = 5\times10^{-4}\). |
| **Heat source** | None (\(q=0\)). |

---

## 2 · Strong form

Find \(T(\mathbf x,t)\) such that

\[
\rho\,c_p\;\frac{\partial T}{\partial t} - \kappa\nabla^{2}T = 0
\quad \text{in}\; \Omega=[0,L]^2,\; t>0
\]

with **periodic boundaries**

\[
T(0,y,t)=T(L,y,t), \qquad 
T(x,0,t)=T(x,L,t),
\]

and initial condition

\[
T(\mathbf x,0)=T_\text{init}(\mathbf x).
\]

Because we prescribe periodicity, there is no Dirichlet (\(\Gamma_D\)) or Neumann (\(\Gamma_N\)) boundary—the domain wraps onto itself.

*You can adjust \(L,\;\alpha,\;\Delta x,\;\Delta t\) via command‑line flags when running the example.*
---

## 3 · Get the demo codes from  [heat2d](files/heat2d.zip)

This folder contains heat2D.c file, a makefile, a post processing file named post2.py and a batch script for running in the cluster

---

## 4 · Quick insights into the code 

### 4.1 Residual assemly

The **Residual()** callback implements the weak form after multiplying the PDE by a test function, integrating by parts, and inserting the B‑spline shape functions.  
For heat conduction the residual per basis function \(N_a\) is  

\[
R_a = -\alpha \, \nabla N_a \!\cdot\! \nabla T ,
\]

coded as:


```c
PetscErrorCode Residual(IGAPoint pnt,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *Re,void *ctx)
{

	AppCtx *user = (AppCtx *)ctx;
  PetscScalar sol; 
  PetscScalar grad[2];
  IGAPointFormValue(pnt,U,&sol);
  IGAPointFormGrad (pnt,U,&grad[0]);
  PetscReal alpha = user->alpha;	
 

  const PetscReal *N0,(*N1)[2];
  IGAPointGetShapeFuns(pnt,0,(const PetscReal**)&N0);
  IGAPointGetShapeFuns(pnt,1,(const PetscReal**)&N1);

  PetscScalar (*Ra)[1] = (PetscScalar (*)[1])Re;
  PetscInt a,nen = pnt->nen;
  for (a=0; a<nen; a++) { Ra[a][0] = -alpha*(N1[a][0]*grad[0] + N1[a][1]*grad[1]); }	

  return 0;
}
```

### 4.2  Initial condition

Modify the hotspot profile in FormInitialCondition():
```c
    T = 100.0*(1-(dist/(0.5*user->Lx*1.414))*(dist/(0.5*user->Lx*1.414)));
```

### 4.3  Output Monitor

OutputMonitor controls the frequency of writing the .dat file which contains the control variables. It also computes the lumped mass matrix vector which is substitued to ts->vec_lump.


### 4.4  Geometry an boundary setup

```c
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetDof(iga,1);CHKERRQ(ierr);

  // setting boundary conditions
  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axis0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,user.Nx,0.0,user.Lx,k);CHKERRQ(ierr);
  IGAAxis axis1;
  ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
  ierr = IGAAxisCopy(axis0,axis1);CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axis1,PETSC_TRUE);CHKERRQ(ierr);

  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
```


### 4.5  Time-stepping driver (explicit Euler)

```c
  /* RHS only */
  ierr = IGASetFormRHSFunction(iga, Residual, &user);CHKERRQ(ierr);

  /* lumped mass vector placeholder */
  ierr = IGACreateVec(iga,&user.lump);
  ierr = VecSet(user.lump,1.0);

  /* create explicit TS object */
  TS ts;
  ierr = IGACreateTS3(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts, TS_LINEAR);
  ierr = TSSetMaxTime(ts,user.total_time);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSEULER);CHKERRQ(ierr);


  /* output monitor */ 
  if (output)  {ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);}
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  /* initial condition and solve */
  Vec C;
  ierr = TSGetSolution(ts,&C);CHKERRQ(ierr);
  ierr = FormInitialCondition(iga,C,&user);CHKERRQ(ierr);
  ierr = TSSolve(ts,C);CHKERRQ(ierr);

```
  

---

## 5 · Build instructions

1. Copy the files into the **`demo/heat/`** folder of your PetIGA clone (or any folder in `$PETIGA_DIR/demo`):

   ```bash
        cp r heat2d  $PETIGA_DIR/demo/
        cd $PETIGA_DIR/demo/heat2d
   ```

2. Compile:

    ```bash
        make heat2D
    ```

## 6 . Running the solver

Edit `run_heat2d.sh` to suit your cluster queue (Slurm, PBS, etc)


**Example run on 128 MPI ranks**

```bash
    #!/bin/bash

    #SBATCH -t 04:00:00
    #SBATCH -n 128
    #SBATCH -o "%x.o%j"
    #SBATCH -e "%x.e%j"
    #SBATCH --job-name="heat2D"
    #SBATCH --mem-per-cpu=1G
    #SBATCH -p RM


    # Run the main program
    mpirun -np 128 ./heat2D > "${SLURM_JOB_NAME}.o$id"
```



## 7 · Results & visualisation *(placeholder)*

**Expected output files**

* `heat2d0000.dat`, `heat2d0100.dat`, … (PetIGA binary snapshots)  
* Energy / time monitor printed to **stdout**

```python
from igakit.io  import PetIGA
from igakit.plot import plt

# load snapshot number 5000 (edit as needed)
xyz, T = PetIGA().read("heat2d0500.dat")

plt.contourf(xyz, T)   # filled contour of temperature
plt.colorbar(label="Temperature")
plt.title("2‑D Heat‑Conduction – step 5000")
plt.show()
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