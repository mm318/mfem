#ifndef CMFEM_ZIENKIEWICZ_ZHU_ESTIMATOR_H
#define CMFEM_ZIENKIEWICZ_ZHU_ESTIMATOR_H

#include "common.h"

CMFEM_STORAGE(ZienkiewiczZhuEstimator, 144);

CMFEM_BEGIN_EXTERN_C

CMFEM_ZienkiewiczZhuEstimator *
CMFEM_ZienkiewiczZhuEstimator_NewDiffusionIntegratorGridFunctionFESpace(
   CMFEM_DiffusionIntegrator *integrator,
   CMFEM_GridFunction *solution,
   CMFEM_FiniteElementSpace *flux_fespace);
CMFEM_ZienkiewiczZhuEstimator *
CMFEM_ZienkiewiczZhuEstimator_NewElasticityIntegratorGridFunctionFESpace(
   CMFEM_ElasticityIntegrator *integrator,
   CMFEM_GridFunction *solution,
   CMFEM_FiniteElementSpace *flux_fespace);
void CMFEM_ZienkiewiczZhuEstimator_Delete(
   CMFEM_ZienkiewiczZhuEstimator *estimator);
void CMFEM_ZienkiewiczZhuEstimator_SetFluxAveraging(
   CMFEM_ZienkiewiczZhuEstimator *estimator,
   int flux_averaging);
void CMFEM_ZienkiewiczZhuEstimator_SetAnisotropic(
   CMFEM_ZienkiewiczZhuEstimator *estimator,
   int anisotropic);

CMFEM_END_EXTERN_C

#endif
