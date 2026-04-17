#ifndef CMFEM_THRESHOLD_REFINER_H
#define CMFEM_THRESHOLD_REFINER_H

#include "common.h"

CMFEM_STORAGE(ThresholdRefiner, 136);

CMFEM_BEGIN_EXTERN_C

CMFEM_ThresholdRefiner *
CMFEM_ThresholdRefiner_NewZienkiewiczZhuEstimator(
   CMFEM_ZienkiewiczZhuEstimator *estimator);
CMFEM_ThresholdRefiner *
CMFEM_ThresholdRefiner_NewLSZienkiewiczZhuEstimator(
   CMFEM_LSZienkiewiczZhuEstimator *estimator);
CMFEM_ThresholdRefiner *
CMFEM_ThresholdRefiner_NewKellyErrorEstimator(
   CMFEM_KellyErrorEstimator *estimator);
void CMFEM_ThresholdRefiner_Delete(CMFEM_ThresholdRefiner *refiner);
void CMFEM_ThresholdRefiner_SetTotalErrorFraction(CMFEM_ThresholdRefiner *refiner,
                                                  double fraction);
void CMFEM_ThresholdRefiner_SetLocalErrorGoal(CMFEM_ThresholdRefiner *refiner,
                                              double goal);
void CMFEM_ThresholdRefiner_PreferConformingRefinement(
   CMFEM_ThresholdRefiner *refiner);
void CMFEM_ThresholdRefiner_SetNCLimit(CMFEM_ThresholdRefiner *refiner,
                                       int nc_limit);
void CMFEM_ThresholdRefiner_Reset(CMFEM_ThresholdRefiner *refiner);
int CMFEM_ThresholdRefiner_Apply(CMFEM_ThresholdRefiner *refiner,
                                 CMFEM_Mesh *mesh);
int CMFEM_ThresholdRefiner_Stop(const CMFEM_ThresholdRefiner *refiner);

CMFEM_END_EXTERN_C

#endif
