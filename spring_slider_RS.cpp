
#include <petscdm.h>
#include <petsc/private/tsimpl.h>
#include <fstream>
#include <iomanip> 

#include "DieterichRuinaAgeing.h"

DieterichRuinaAgeing alwa;
std::ofstream out_file;


static PetscErrorCode RHSFunction_spring_slider(TS ts, PetscReal t, Vec U, Vec F, void *ctx)
{
  PetscScalar       *f;
  const PetscScalar *u;
  PetscScalar       D, psi;

  
  DieterichRuinaAgeing* alwa = static_cast<DieterichRuinaAgeing*>(ctx);
  double tau;
  
  double V;



  
  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));
  PetscCall(VecGetArray(F, &f));

  D = u[0];
  psi = u[1];
  tau = alwa->k * ((alwa->Vp * t + alwa->yield_point_init) - D);
  V = alwa->slip_rate(tau, psi);

  f[0] = V;
  f[1] = alwa->state_rhs(V, psi);

  PetscCall(VecRestoreArrayRead(U, &u));
  PetscCall(VecRestoreArray(F, &f));
  PetscFunctionReturn(0);
}


PetscErrorCode ts_soln_view(TS ts)
{
  Vec U;
  const PetscScalar *u;
  PetscInt step;
  PetscScalar time;
  PetscFunctionBeginUser;
  PetscCall(TSGetStepNumber(ts, &step));
  PetscCall(TSGetTime(ts, &time));
  PetscCall(TSGetSolution(ts, &U));
  PetscCall(VecGetArrayRead(U, &u));
  PetscScalar D = u[0];
  PetscScalar psi = u[1];
  double tau = alwa.k * (alwa.Vp * time - D);
  double V = alwa.slip_rate(tau, psi);
  out_file << std::scientific << std::setprecision(4) << (double)time << "," << (double)D << "," << (double)psi << "," << V << "," << tau << std::endl;
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscFunctionReturn(0);
}


int main(int argc, char **argv)
{
  if (argc != 14) {
    std::cout << "Usage: " << argv[0] << " V0 f0 a b eta L sn Vinit Vp k yield_point_init final_time" << std::endl;
    return 1;
  }
  alwa.V0 = std::stof(argv[1]);
  alwa.f0 = std::stof(argv[2]);
  alwa.a = std::stof(argv[3]);
  alwa.b = std::stof(argv[4]);
  alwa.eta = std::stof(argv[5]);
  alwa.L = std::stof(argv[6]);
  alwa.sn = std::stof(argv[7]);
  alwa.Vinit = std::stof(argv[8]);
  alwa.Vp = std::stof(argv[9]);
  alwa.k = std::stof(argv[10]);
  alwa.yield_point_init = std::stof(argv[11]);


  double finle_time = std::stof(argv[12]);
  double tau_init = alwa.k * alwa.yield_point_init;
  double psi_init = alwa.psi_init(tau_init);

  std::string out_file_name = argv[13];
  


  std::cout << "Parameters:" << std::endl;
  std::cout << "V0 = " << alwa.V0 << std::endl;
  std::cout << "f0 = " << alwa.f0 << std::endl;
  std::cout << "a = " << alwa.a << std::endl;
  std::cout << "b = " << alwa.b << std::endl;
  std::cout << "eta = " << alwa.eta << std::endl;
  std::cout << "L = " << alwa.L << std::endl;
  std::cout << "sn = " << alwa.sn << std::endl;
  std::cout << "Vinit = " << alwa.Vinit << std::endl;
  std::cout << "Vp = " << alwa.Vp << std::endl;
  std::cout << "k = " << alwa.k << std::endl;
  std::cout << "Initial state:" << psi_init <<std::endl;
  std::cout << "final_time = " << finle_time << std::endl;

  out_file.open(out_file_name);
  if (!out_file) {
      std::cerr << "Error opening out.txt" << std::endl;
      return 1;
      }
  out_file << "t,D,psi,V,tau" << std::endl;

  TS           ts; /* ODE integrator */
  Vec          U;  /* solution will be stored here */
  PetscMPIInt  size;
  PetscInt     n = 2;
  PetscScalar *u;
  TSAdapt      adapt;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetType(ts, TSRK));
  PetscCall(TSSetProblemType(ts, TS_NONLINEAR));

  PetscCall(TSSetRHSFunction(ts, NULL, RHSFunction_spring_slider, static_cast<void*>(&alwa)));

  /* initial condition */
  PetscCall(VecCreate(PETSC_COMM_WORLD, &U));
  PetscCall(VecSetSizes(U, n, PETSC_DETERMINE));
  PetscCall(VecSetUp(U));
  PetscCall(VecGetArray(U, &u));
  u[0] = 0.0;
  u[1] = psi_init;
  PetscCall(VecRestoreArray(U, &u));

  PetscCall(TSSetSolution(ts, U));

  PetscCall(TSSetMaxTime(ts, finle_time));
  PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER));
  PetscCall(TSSetTimeStep(ts, 1.0e-10));
  /* The adaptive time step controller could take very
     large timesteps. An upper limit is enforced here to avoid this. */
  PetscCall(TSGetAdapt(ts, &adapt));
  PetscCall(TSAdaptSetStepLimits(adapt, 0.0, 2000000.0));

  PetscCall(TSSetPostStep(ts, ts_soln_view));
  
  PetscCall(TSSetFromOptions(ts));

  PetscCall(TSSetUp(ts));

  PetscCall(TSSolve(ts, U));
  
  PetscCall(VecDestroy(&U));
  PetscCall(TSDestroy(&ts));

  PetscCall(PetscFinalize());

  out_file.close();
  return 0;
}
