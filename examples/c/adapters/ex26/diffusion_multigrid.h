#ifndef CMFEM_DIFFUSION_MULTIGRID_H
#define CMFEM_DIFFUSION_MULTIGRID_H

#include "../../../../c_api/common.h"

CMFEM_STORAGE(DiffusionMultigrid, 40);

CMFEM_BEGIN_EXTERN_C

CMFEM_DiffusionMultigrid *CMFEM_DiffusionMultigrid_NewMeshGrOrAi(
   CMFEM_Mesh *mesh,
   int geometric_refinements,
   int order_refinements,
   const CMFEM_ArrayInt *ess_bdr);
void CMFEM_DiffusionMultigrid_Delete(CMFEM_DiffusionMultigrid *multigrid);
CMFEM_FiniteElementSpace *CMFEM_DiffusionMultigrid_GetFinestFESpace(
   CMFEM_DiffusionMultigrid *multigrid);
CMFEM_Mesh *CMFEM_DiffusionMultigrid_GetFinestMesh(
   CMFEM_DiffusionMultigrid *multigrid);
void CMFEM_DiffusionMultigrid_SetCycleTypeV(
   CMFEM_DiffusionMultigrid *multigrid,
   int pre_smoothing,
   int post_smoothing);
void CMFEM_DiffusionMultigrid_FormFineLinearSystemOp(
   CMFEM_DiffusionMultigrid *multigrid,
   CMFEM_GridFunction *x,
   CMFEM_LinearForm *b,
   CMFEM_OperatorPtr *A,
   CMFEM_Vector *X,
   CMFEM_Vector *B);
void CMFEM_DiffusionMultigrid_PCGOp(
   CMFEM_DiffusionMultigrid *multigrid,
   const CMFEM_OperatorPtr *A,
   const CMFEM_Vector *B,
   CMFEM_Vector *X,
   int print_iter,
   int max_iter,
   double rtol,
   double atol);
void CMFEM_DiffusionMultigrid_RecoverFineFEMSolution(
   CMFEM_DiffusionMultigrid *multigrid,
   const CMFEM_Vector *X,
   const CMFEM_LinearForm *b,
   CMFEM_GridFunction *x);

CMFEM_END_EXTERN_C

#endif
