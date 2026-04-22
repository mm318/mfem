#ifndef CMFEM_ELASTICITY_INTEGRATOR_H
#define CMFEM_ELASTICITY_INTEGRATOR_H

#include "../common.h"

CMFEM_STORAGE(ElasticityIntegrator, 368);

CMFEM_BEGIN_EXTERN_C

CMFEM_ElasticityIntegrator *CMFEM_ElasticityIntegrator_NewPwcPwc(
   const CMFEM_PWConstCoefficient *lambda,
   const CMFEM_PWConstCoefficient *mu);
void CMFEM_ElasticityIntegrator_Delete(CMFEM_ElasticityIntegrator *integrator);

CMFEM_END_EXTERN_C

#endif
