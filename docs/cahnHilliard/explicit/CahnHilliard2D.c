#include "petiga.h"
#include <petsc/private/tsimpl.h>


typedef struct {
  IGA  iga;
  PetscReal theta,alpha;
  PetscReal cbar;
  PetscReal Eprev;
  PetscInt  step;
  Vec lump;
} AppCtx;

void Mobility(AppCtx *user,PetscReal c,PetscReal *M,PetscReal *dM,PetscReal *d2M)
{
  if (M)   *M   = c*(1-c);
  if (dM)  *dM  = 1-2*c;
  if (d2M) *d2M = -2;
}

void ChemicalPotential(AppCtx *user,PetscReal c,PetscReal *mu,PetscReal *dmu,PetscReal *d2mu)
{
  PetscReal theta  = user->theta;
  PetscReal alpha  = user->alpha;
  if (mu) {
    (*mu)  = 0.5/theta*log(c/(1-c)) + 1 - 2*c;
    (*mu) *= 3*alpha;
  }
  if (dmu) {
    (*dmu)  = 0.5/theta*1/(c*(1-c)) - 2;
    (*dmu) *= 3*alpha;
  }
  if (d2mu) {
    (*d2mu)  = 0.5/theta*(2*c-1)/(c*c*(1-c)*(1-c));
    (*d2mu) *= 3*alpha;
  }
}


PetscReal GinzburgLandauFreeEnergy(PetscReal c,PetscReal cx,PetscReal cy,AppCtx *user)
{
  PetscReal theta = user->theta;
  PetscReal alpha = user->alpha;
  PetscReal E = c*log(c) + (1-c)*log(1-c) + 2*theta*c*(1-c) + theta/(3*alpha)*(cx*cx+cy*cy);
  return E;
}

PetscErrorCode Stats(IGAPoint p,const PetscScalar *U,PetscInt n,PetscScalar *S,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;
  PetscFunctionBegin;

  PetscScalar c,c1[3];
  IGAPointFormValue(p,U,&c);
  IGAPointFormGrad(p,U,&c1[0]);
  PetscReal diff = c - user->cbar;

  S[0] = GinzburgLandauFreeEnergy(c,c1[0],c1[1],user); // Free energy
  S[1] = diff*diff;                                    // Second moment
  S[2] = S[1]*diff;                                    // Third moment

  PetscFunctionReturn(0);
}

PetscErrorCode StatsMonitor(TS ts,PetscInt step,PetscReal t,Vec U,void *mctx)
{
  AppCtx         *user = (AppCtx *)mctx;
  IGA            iga;
  PetscReal      dt;
  PetscScalar    stats[3] = {0,0,0};
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)ts,"IGA",(PetscObject*)&iga);CHKERRQ(ierr);
  PetscValidHeaderSpecific(iga,IGA_CLASSID,0);

  ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
  ierr = IGAComputeScalar(iga,U,3,&stats[0],Stats,mctx);CHKERRQ(ierr);

  if (step == 0) {ierr = PetscPrintf(PETSC_COMM_WORLD,"#Step    Time       dt           Free Energy            Second moment          Third moment\n");CHKERRQ(ierr);}
  if(step%1000 == 0) {ierr = PetscPrintf(PETSC_COMM_WORLD,"%d,%.6e %.6e %.16e %.16e %.16e\n",(int)step,(double)t,(double)dt,(double)stats[0],(double)stats[1],(double)stats[2]);CHKERRQ(ierr);}

  if (step == 0) user->Eprev = PETSC_MAX_REAL;
  if((PetscReal)stats[0] > user->Eprev) {ierr = PetscPrintf(PETSC_COMM_WORLD,"WARNING: Ginzburg-Landau free energy increased!\n");CHKERRQ(ierr);}
  user->Eprev = PetscRealPart(stats[0]);
  PetscFunctionReturn(0);
}

PetscErrorCode Residual(IGAPoint p,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *R,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

  PetscInt nen;
  IGAPointGetSizes(p,0,&nen,0);

  PetscScalar c_t,c;
  //IGAPointFormValue(p,V,&c_t);
  IGAPointFormValue(p,U,&c);

  PetscReal M,dM;
  Mobility(user,c,&M,&dM,NULL);
  PetscReal dmu;
  ChemicalPotential(user,c,NULL,&dmu,NULL);

  PetscScalar c1[2],del2_c;
  IGAPointFormGrad(p,U,&c1[0]);
  IGAPointFormDel2(p,U,&del2_c);
  PetscScalar c_x = c1[0], c_y = c1[1];

  PetscScalar t1 = M*dmu + dM*del2_c;

  const PetscReal (*N0)       = (typeof(N0)) p->shape[0];
  const PetscReal (*N1)[2]    = (typeof(N1)) p->shape[1];
  const PetscReal (*N2)[2][2] = (typeof(N2)) p->shape[2];

  PetscInt a;
  for (a=0; a<nen; a++) {
    PetscReal Na    = N0[a];
    PetscReal Na_x  = N1[a][0];
    PetscReal Na_y  = N1[a][1];
    PetscReal Na_xx = N2[a][0][0];
    PetscReal Na_yy = N2[a][1][1];
    /* ----- */
    PetscScalar Ra  = 0;
    // grad(Na) . ((M*dmu + dM*del2(c))) grad(C)
    Ra -= (Na_x * c_x + Na_y * c_y) * t1;
    // del2(Na) * M * del2(c)
    Ra -= (Na_xx+Na_yy) * M * del2_c;
    /* ----- */
    R[a] = Ra;
  }
  return 0;
}


PetscErrorCode FormInitialCondition(IGA iga,Vec C,AppCtx *user)
{
  MPI_Comm       comm;
  PetscRandom    rctx;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
  ierr = PetscRandomCreate(comm,&rctx);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rctx,user->cbar-0.05,user->cbar+0.05);CHKERRQ(ierr);
  ierr = PetscRandomSeed(rctx);CHKERRQ(ierr);
  ierr = VecSetRandom(C,rctx);CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OutputMonitor(TS ts,PetscInt step,PetscReal t,Vec U,void *mctx)
{
  AppCtx         *user = (AppCtx *)mctx;
  IGA            iga;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)ts,"IGA",(PetscObject*)&iga);CHKERRQ(ierr);
  PetscValidHeaderSpecific(iga,IGA_CLASSID,0);


  if(user->step==0){

    IGAElement        element;
    IGAPoint          point;
    void              *ctx;
    PetscScalar       *J,*K;

    PetscFunctionBegin;
    PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
    IGACheckSetUp(iga,1);

    Mat		matJtemp;
    ierr = IGACreateMat(iga,&matJtemp);
    ierr = MatZeroEntries(matJtemp);CHKERRQ(ierr);

    ierr = IGABeginElement(user->iga,&element);CHKERRQ(ierr);
    while (IGANextElement(user->iga,element)) {
      ierr = IGAElementGetWorkMat(element,&J);CHKERRQ(ierr);
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
        ierr = IGAPointGetWorkMat(point,&K);CHKERRQ(ierr);
        const PetscReal *N0,(*N1)[2];
        IGAPointGetShapeFuns(point,0,(const PetscReal**)&N0);
        IGAPointGetShapeFuns(point,1,(const PetscReal**)&N1);
        PetscInt dof = point->dof,a,b,c,nen = point->nen;
        PetscScalar (*Je)[dof][nen][dof] = (PetscScalar (*)[dof][nen][dof])K;
        for (a=0; a<nen; a++){
          for (b=0; b<nen; b++){
            for (c=0; c<dof; c++){
              Je[a][c][b][c] += N0[a]*N0[b];
            }
          }
        }
        ierr = IGAPointAddMat(point,K,J);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
      ierr = IGAElementFixJacobian(element,J);CHKERRQ(ierr);
      ierr = IGAElementAssembleMat(element,J,matJtemp);CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga,&element);CHKERRQ(ierr);


    ierr = MatAssemblyBegin(matJtemp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (matJtemp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    Vec		    temp;
    ierr = IGACreateVec(iga,&temp);
    ierr = VecSet(temp,1.0);
    ierr = MatMultTranspose(matJtemp,temp,user->lump);

    VecDestroy(&temp);
    MatDestroy(&matJtemp);

    ts->vec_lump = user->lump;

    char filelump[256];
    sprintf(filelump, "./lump.dat"); 
    ierr = IGAWriteVec(user->iga,user->lump,filelump);CHKERRQ(ierr);

    PetscScalar *arrayLC;
    PetscInt  b, NlocLC;
    ierr = VecGetArray(user->lump,&arrayLC);CHKERRQ(ierr);
    ierr = VecGetLocalSize(user->lump,&NlocLC);CHKERRQ(ierr);
    ierr = VecRestoreArray(user->lump,&arrayLC);CHKERRQ(ierr);

  }


  if(user->step%100000 == 0){
    ierr = PetscSNPrintf(filename,sizeof(filename),"./ch2d%d.dat",(int)step);CHKERRQ(ierr);
    ierr = IGAWriteVec(iga,U,filename);CHKERRQ(ierr);
  }
  user->step++;
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {

  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  /* Define simulation specific parameters */
  AppCtx user;
  user.alpha = 3000.0; /* interface thickess parameter */
  user.theta = 1.5;    /* temperature/critical temperature */
  user.cbar  = 0.63;   /* average concentration */
  user.step = 0;

  /* Set discretization options */
  PetscInt  N = 100;
  PetscInt  p = 2;
  PetscInt  k = PETSC_DECIDE;
  char      initial[PETSC_MAX_PATH_LEN] = {0};
  PetscBool output  = PETSC_TRUE;
  PetscBool monitor = PETSC_TRUE;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","CahnHilliard2D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-N","number of elements (along one dimension)",__FILE__,N,&N,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-p","polynomial order",__FILE__,p,&p,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-k","global continuity order",__FILE__,k,&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-initial","Load initial solution from file",__FILE__,initial,initial,sizeof(initial),NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-output","Enable output files",__FILE__,output,&output,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-monitor","Compute and show statistics of solution",__FILE__,monitor,&monitor,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-cbar","Initial average concentration",__FILE__,user.cbar,&user.cbar,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-alpha","Interface thickess parameter",__FILE__,user.alpha,&user.alpha,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-theta","Ratio temperature/critical temperature",__FILE__,user.theta,&user.theta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (k == PETSC_DECIDE) k = p-1;

  if (p < 2 || k < 1) /* Problem requires a p>=2 C^1 basis */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,
            "Problem requires minimum of p = 2 and k = 1");
  if (p <= k)         /* Check k < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,
            "Discretization inconsistent: polynomial order must be greater than degree of continuity");

  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetDof(iga,1);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axis0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,N,0.0,1.0,k);CHKERRQ(ierr);
  IGAAxis axis1;
  ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
  ierr = IGAAxisCopy(axis0,axis1);CHKERRQ(ierr);

  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  ierr = IGAWrite(iga,"igaphase.dat");CHKERRQ(ierr);
  user.iga = iga;


  ierr = IGASetFormRHSFunction(iga,Residual,&user);CHKERRQ(ierr);

  ierr = IGACreateVec(iga,&user.lump);
  ierr = VecSet(user.lump,1.0);


  TS ts;
  ierr = IGACreateTS3(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts, TS_LINEAR);
  ierr = TSSetMaxTime(ts,1.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,1e-9);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSEULER);CHKERRQ(ierr);

  if (monitor) {ierr = TSMonitorSet(ts,StatsMonitor,&user,NULL);CHKERRQ(ierr);}
  if (output)  {ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);}
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  Vec C;
  ierr = TSGetSolution(ts,&C);CHKERRQ(ierr);
  if (initial[0] != 0) { /* initial condition from datafile */
    ierr = IGAReadVec(iga,C,initial);CHKERRQ(ierr);
  } else {                /* initial condition is random */
    ierr = FormInitialCondition(iga,C,&user);CHKERRQ(ierr);
  }



  ierr = TSSolve(ts,C);CHKERRQ(ierr);

  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
