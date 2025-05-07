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

*(You can adjust \(L,\;\alpha,\;\Delta x,\;\Delta t\) via command‑line flags when running the example.)*
---

## 3 · Get the demo codes

| Variant | Download link |
|---------|---------------|
| **Implicit PetIGA solver** | [heat2d\_implicit.c](files/heat2d_implicit.c) |
| **Explicit PetIGA solver** | [heat2d\_explicit.c](files/heat2d_explicit.c) |
| **Sample `run.sh` script** | [run\_heat2d.sh](files/run_heat2d.sh) |

*(Replace the links with actual files in `docs/files/`.)*

---

## 4 · Build instructions

1. Copy the files into the **`demo/heat/`** folder of your PetIGA clone (or any folder in `$PETIGA_DIR/demo`):

   ```bash
        cp heat2d_*.c  $PETIGA_DIR/demo/heat/
        cp run_heat2d.sh $PETIGA_DIR/demo/heat/
        cd $PETIGA_DIR/demo/heat
   ```

2. Compile:

    ```bash
        # Pick one; or build both
        make heat2d_implicit
        make heat2d_explicit
    ```

## 5 . Running the solver

Edit `run_heat2d.sh` to uit your cluster queue (Slurm, PBS, etc)

Key options you migh change:

| Option        | Meaning                                 | Default |
|---------------|-----------------------------------------|---------|
| `-N`          | number of elements per side             | 32      |
| `-p`          | polynomial degree                       | 2       |
| `-ts_type`    | time‑stepping scheme (`euler`, `rk3`, `theta`, …) | `euler` (explicit build) |
| `-ts_dt`      | time step                               | 1e‑4    |
| `-ts_max_time`| final simulation time                   | 0.5     |

**Example run on 4 MPI ranks**

```bash
mpirun -np 4 ./heat2d_explicit \
       -N 64          \
       -ts_dt 2e-4    \
       -ts_max_time 1.0
# add extra flags if needed, e.g. -iga_view or solver tolerances
```



## 6 · Results & visualisation *(placeholder)*

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