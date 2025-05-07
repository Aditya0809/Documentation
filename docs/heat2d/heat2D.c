#include "petiga.h"
#include <petsc/private/tsimpl.h>

#define SQ(x) ((x)*(x))


// USER DEFINED STRUCTURE
typedef struct {
  IGA  iga;
  PetscInt Nx,Ny,Lx,Ly, nprintdat, nprintu;
  PetscInt p,k,count;
  PetscReal alpha,dt,total_time;
  Vec lump;
} AppCtx;


typedef struct {
 PetscScalar T;
} Field;






// FUNCTION TO DEFINE THE RHS FUNCTION VECTOR
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




/* FUNCTION TO DEFINE THE INTIAL CONDITION */
PetscErrorCode FormInitialCondition(IGA iga,Vec U,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
 
 
  
  DM da;
  ierr = IGACreateNodeDM(iga,1,&da);CHKERRQ(ierr);
  Field **u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  PetscReal pi = 3.14159265358979323846;


  PetscInt i,j,k;
  PetscInt xs = info.xs, ys = info.ys, xm = info.xm, ym = info.ym;
  for(i = info.xs; i < info.xs + info.xm; i++)
  {
    for(j = info.ys; j < info.ys + info.ym; j++)
    { 
      PetscReal x = (PetscReal)i*user->Lx/(user->Nx-1); 
      PetscReal y = (PetscReal)j*user->Lx/(user->Nx-1); 
      PetscReal dist,T;
      dist = sqrt((x - 0.5*user->Lx)*(x - 0.5*user->Lx) + (y - 0.5*user->Lx)*(y - 0.5*user->Lx));
      T = 100.0*(1-(dist/(0.5*user->Lx*1.414))*(dist/(0.5*user->Lx*1.414)));
      u[j][i].T = T;
    }
  }
  
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}






PetscErrorCode OutputMonitor(TS ts,PetscInt step,PetscReal t,Vec U,void *mctx)
{
  AppCtx         *user = (AppCtx *)mctx;
  IGA            iga;
  char           filename[256];
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)ts,"IGA",(PetscObject*)&iga);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject)U,"IGA",(PetscObject*)&iga);CHKERRQ(ierr);
  PetscValidHeaderSpecific(iga,IGA_CLASSID,0);

  // creating the  lumped mass matrix vector just once
  if(step==0){

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

  PetscInt fac = 2.0/user->dt;

  user->count = 0;

  if(step%fac==0)
  { 
    // PRINTING OUTPUT FILES
	  sprintf(filename,"./heat%d.dat",step); 
	  ierr = IGAWriteVec(iga,U,filename);CHKERRQ(ierr);
  }


  PetscFunctionReturn(0);
}








int main(int argc, char *argv[]) {

  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

 /* Define simulation specific parameters */
  AppCtx user; 
  user.alpha = 5.0; // Diffusion coefficient
  user.dt = 0.0005; // delta t
  user.Nx = 100; // number of elements in x direction = y direction
  user.Lx = 50.0; // length of the domain
  user.count = 0; 
  user.total_time = 50.0; // total time of simulation

  /* Set discretization options */

  PetscInt  N;
  PetscReal L, total_time, dt, alpha;
  PetscInt  p = 2; // polynomial order
  PetscInt  k = PETSC_DECIDE;
  char      initial[PETSC_MAX_PATH_LEN] = {0};
  PetscBool output  = PETSC_TRUE;
  PetscBool monitor = PETSC_FALSE;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","Heat conduction 2D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-N","number of elements (along one dimension)",__FILE__,N,&N,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-p","polynomial order",__FILE__,p,&p,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-k","global continuity order",__FILE__,k,&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-L","length of the domain",__FILE__,L,&L,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-total_time","total time of simulation",__FILE__,total_time,&total_time,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-dt","time step",__FILE__,dt,&dt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-alpha","Diffusion coefficient",__FILE__,alpha,&alpha,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-initial","Load initial solution from file",__FILE__,initial,initial,sizeof(initial),NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-output","Enable output files",__FILE__,output,&output,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-monitor","Compute and show statistics of solution",__FILE__,monitor,&monitor,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (k == PETSC_DECIDE) k = p-1;

  user.Nx = N;
  user.Lx = L;
  user.total_time = total_time;
  user.dt = dt;
  user.alpha = alpha;

  if (p < 2 || k < 1) /* Problem requires a p>=2 C^1 basis */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,
            "Problem requires minimum of p = 2 and k = 1");
  if (p <= k)         /* Check k < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,
            "Discretization inconsistent: polynomial order must be greater than degree of continuity");

  // IGA for original system
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
  ierr = IGAWrite(iga,"igaphase.dat");CHKERRQ(ierr);
  user.iga = iga;

  PetscInt dim;
  ierr = IGAGetDim(iga,&dim);CHKERRQ(ierr);

  // Set up the rhs function
  ierr = IGASetFormRHSFunction(iga, Residual, &user);CHKERRQ(ierr);

  ierr = IGACreateVec(iga,&user.lump);
  ierr = VecSet(user.lump,1.0);

  TS ts;
  ierr = IGACreateTS3(iga,&ts);CHKERRQ(ierr);
 
  ierr = TSSetProblemType(ts, TS_LINEAR);
  ierr = TSSetMaxTime(ts,user.total_time);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSEULER);CHKERRQ(ierr);



  if (output)  {ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);}
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  Vec C;
  ierr = TSGetSolution(ts,&C);CHKERRQ(ierr);
  ierr = FormInitialCondition(iga,C,&user);CHKERRQ(ierr);
  ierr = TSSolve(ts,C);CHKERRQ(ierr);
  


  ierr = VecDestroy(&C);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
