#ifndef CMFEM_KELLY_ERROR_ESTIMATOR_H
#define CMFEM_KELLY_ERROR_ESTIMATOR_H

#include "common.h"

CMFEM_STORAGE(KellyErrorEstimator, 224);

CMFEM_BEGIN_EXTERN_C

CMFEM_KellyErrorEstimator *
CMFEM_KellyErrorEstimator_NewDiGfFes(
   CMFEM_DiffusionIntegrator *integrator,
   CMFEM_GridFunction *solution,
   CMFEM_FiniteElementSpace *flux_fespace);
void CMFEM_KellyErrorEstimator_Delete(CMFEM_KellyErrorEstimator *estimator);
void CMFEM_KellyErrorEstimator_Reset(CMFEM_KellyErrorEstimator *estimator);

CMFEM_END_EXTERN_C

#endif
