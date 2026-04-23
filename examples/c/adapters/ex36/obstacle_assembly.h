#ifndef CMFEM_OBSTACLE_ASSEMBLY_H
#define CMFEM_OBSTACLE_ASSEMBLY_H

#include "../../../../c_api/common.h"

CMFEM_BEGIN_EXTERN_C

void CMFEM_ProjectObstacleLogGfFc(CMFEM_GridFunction *grid_function,
                                  const CMFEM_GridFunction *state,
                                  const CMFEM_FunctionCoefficient *obstacle,
                                  double min_val);
void CMFEM_ProjectObstacleExpGfFc(CMFEM_GridFunction *grid_function,
                                  const CMFEM_GridFunction *state,
                                  const CMFEM_FunctionCoefficient *obstacle,
                                  double min_val,
                                  double max_val);
void CMFEM_AssembleObstacleLfGfFc(CMFEM_FiniteElementSpace *fespace,
                                  const CMFEM_GridFunction *state,
                                  const CMFEM_FunctionCoefficient *obstacle,
                                  CMFEM_Vector *out,
                                  double min_val,
                                  double max_val);
CMFEM_SparseMatrix *CMFEM_AssembleObstacleHessianSmGf(
   CMFEM_FiniteElementSpace *fespace,
   const CMFEM_GridFunction *state,
   int order,
   double spectrum_shift,
   double min_val,
   double max_val);

CMFEM_END_EXTERN_C

#endif
