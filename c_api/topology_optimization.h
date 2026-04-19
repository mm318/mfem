#ifndef CMFEM_TOPOLOGY_OPTIMIZATION_H
#define CMFEM_TOPOLOGY_OPTIMIZATION_H

#include "common.h"

CMFEM_BEGIN_EXTERN_C

CMFEM_Mesh *CMFEM_Ex37_NewBeamMesh(void);

void *CMFEM_Ex37FilterSolver_NewMeshOrderEpsAi(CMFEM_Mesh *mesh,
                                               int order,
                                               double epsilon_sq,
                                               const CMFEM_ArrayInt *ess_bdr);
void CMFEM_Ex37FilterSolver_Delete(void *solver);
void CMFEM_Ex37FilterSolver_SolveSigmoidGf(void *solver,
                                           const CMFEM_GridFunction *psi,
                                           CMFEM_GridFunction *rho_filter);
void CMFEM_Ex37FilterSolver_SolveSedGfGf(void *solver,
                                         double lambda,
                                         double mu,
                                         const CMFEM_GridFunction *u,
                                         const CMFEM_GridFunction *rho_filter,
                                         double rho_min,
                                         CMFEM_GridFunction *w_filter);

void *CMFEM_Ex37ElasticitySolver_NewMeshOrderAiCircle(
   CMFEM_Mesh *mesh,
   int order,
   const CMFEM_ArrayInt *ess_bdr,
   double radius,
   const double *center_xy,
   const double *force_xy);
void CMFEM_Ex37ElasticitySolver_Delete(void *solver);
void CMFEM_Ex37ElasticitySolver_SolveGf(void *solver,
                                        double lambda,
                                        double mu,
                                        const CMFEM_GridFunction *rho_filter,
                                        double rho_min,
                                        CMFEM_GridFunction *u);
double CMFEM_Ex37ElasticitySolver_ComplianceGf(void *solver,
                                               const CMFEM_GridFunction *u);

double CMFEM_Ex37_ProjectSigmoidVolumeGf(CMFEM_GridFunction *psi,
                                         double target_volume,
                                         double tol,
                                         int max_its);
double CMFEM_Ex37_SigmoidDiffL1GfGf(const CMFEM_GridFunction *psi,
                                    const CMFEM_GridFunction *psi_old);
void CMFEM_Ex37_ProjectSigmoidGf(const CMFEM_GridFunction *psi,
                                 CMFEM_GridFunction *rho_gf);
void CMFEM_Ex37_ProjectSimpGf(const CMFEM_GridFunction *rho_filter,
                              double rho_min,
                              CMFEM_GridFunction *rho_gf);

CMFEM_END_EXTERN_C

#endif
