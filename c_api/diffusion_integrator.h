#ifndef CMFEM_DIFFUSION_INTEGRATOR_H
#define CMFEM_DIFFUSION_INTEGRATOR_H

#include "common.h"

CMFEM_STORAGE(DiffusionIntegrator, 1008);

CMFEM_BEGIN_EXTERN_C

CMFEM_DiffusionIntegrator *CMFEM_DiffusionIntegrator_NewCc(
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_DiffusionIntegrator_Delete(CMFEM_DiffusionIntegrator *integrator);

CMFEM_END_EXTERN_C

#endif
