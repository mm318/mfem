#ifndef CMFEM_LS_ZIENKIEWICZ_ZHU_ESTIMATOR_H
#define CMFEM_LS_ZIENKIEWICZ_ZHU_ESTIMATOR_H

#include "common.h"

CMFEM_STORAGE(LSZienkiewiczZhuEstimator, 104);

CMFEM_BEGIN_EXTERN_C

CMFEM_LSZienkiewiczZhuEstimator *
CMFEM_LSZienkiewiczZhuEstimator_NewDiffusionIntegratorGridFunction(
   CMFEM_DiffusionIntegrator *integrator,
   CMFEM_GridFunction *solution);
void CMFEM_LSZienkiewiczZhuEstimator_Delete(
   CMFEM_LSZienkiewiczZhuEstimator *estimator);
void CMFEM_LSZienkiewiczZhuEstimator_SetTichonovRegularization(
   CMFEM_LSZienkiewiczZhuEstimator *estimator,
   double coefficient);

CMFEM_END_EXTERN_C

#endif
