# Tutorial II · 2‑D Cahn Hilliard Equation
*(explicit **vs** implicit PetIGA implementations)*


In this example we focus on:

1. The mathematical model (strong form)  
2. How the **Residual** and **Tangent** functions differ between explicit and implicit time stepping  
3. Minimal IGA setup common to both variants  
4. Place-holders for running, post-processing, and comparing results


## 1 · Problem description 

* 2-D periodic square Ω = [0,1]²  
* Order-parameter \(c(\mathbf x,t)\) representing concentration (0 ≤ c ≤ 1)  
* Temperature ratio θ = 1.5, interface-parameter α = 3000  
* Mobility \(M(c)=c(1-c)\)  
* Random initial condition with mean \( \bar c = 0.63\)

*(Feel free to expand.)*

---



## 2 · Strong form

The Cahn–Hilliard system couples a **conserved** order parameter \(c\) with its chemical potential \(\mu\):

\[
\begin{aligned}
\frac{\partial c}{\partial t} &= \nabla\!\cdot\!\bigl(M(c)\,\nabla\mu \bigr), 
\\[4pt]
\mu &= \frac{\partial\Psi}{\partial c} \;-\; \alpha\,\theta\,\nabla^{2}c,
\end{aligned}
\]

where  

\[
\Psi(c) = \tfrac12 \theta\,\bigl[c\ln c + (1-c)\ln(1-c)\bigr] + \theta\,c(1-c)
\]

yields  

\[
\frac{\partial\Psi}{\partial c}
\;=\;
\tfrac12\theta \ln\!\left(\frac{c}{1-c}\right)+\theta\,(1-2c).
\]

Boundary conditions are **periodic**; source term \(q=0\).


## 3 · Code structure overview

| File | Purpose |
|------|---------|
| `ch2d_explicit.c` | explicit Euler (`TSEULER`) – **Residual only** |
| `ch2d_implicit.c` | implicit (`TSALPHA`, or `TSTheta`) – **Residual + Tangent** |
| `run_ch2d.sh` | batch script template |

> We show only the parts that differ between the two builds—everything else (domain, axis settings, initial condition, monitors) is shared.

---


### 3.1 Explicit variant – Residual only

```c
/* ---------------- Residual for explicit solver ---------------- */
PetscErrorCode Residual(IGAPoint p,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *R,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

  /* c,  ∇c,  ∇²c  */
  PetscScalar  c, gradc[2], del2c;
  IGAPointFormValue (p,U,&c);
  IGAPointFormGrad  (p,U,&gradc[0]);
  IGAPointFormDel2  (p,U,&del2c);

  /* mobility and d(mu)/dc */
  PetscReal M,dM;
  Mobility(user,c,&M,&dM,NULL);

  PetscReal dmu;
  ChemicalPotential(user,c,NULL,&dmu,NULL);

  /* weak form */
  const PetscReal (*N0)    = (typeof(N0)) p->shape[0];
  const PetscReal (*N1)[2] = (typeof(N1)) p->shape[1];
  const PetscReal (*N2)[2][2] = (typeof(N2)) p->shape[2];

  PetscInt a, nen = p->nen;
  for (a=0; a<nen; a++) {
    PetscReal Na_x  = N1[a][0];
    PetscReal Na_y  = N1[a][1];
    PetscReal Na_xx = N2[a][0][0];
    PetscReal Na_yy = N2[a][1][1];

    PetscScalar Ra  = 0.0;
    PetscScalar t1  = M*dmu + dM*del2c;

    Ra -= (Na_x*gradc[0] + Na_y*gradc[1]) * t1;     /* convective part   */
    Ra -= (Na_xx+Na_yy) * M * del2c;                /* diffusion part    */
    R[a] = Ra;
  }
  return 0;
}
```


No Jacobian is registered; the driver simply calls

```c
    IGASetFormRHSFunction(iga, Residual, &user);
    IGACreateTS3(iga,&ts);      /* attaches RHS only          */
    TSSetType(ts, TSEULER);     /* forward Euler time stepper */
```
