#ifndef CMFEM_HYPERELASTIC_DYNAMICS_H
#define CMFEM_HYPERELASTIC_DYNAMICS_H

#include "common.h"

CMFEM_BEGIN_EXTERN_C

void *CMFEM_Ex10HyperelasticOperator_New(CMFEM_FiniteElementSpace *fespace,
                                         const CMFEM_ArrayInt *ess_bdr,
                                         double viscosity,
                                         double mu,
                                         double bulk_modulus);
void CMFEM_Ex10HyperelasticOperator_Delete(void *oper);
int CMFEM_Ex10HyperelasticOperator_StateSize(const void *oper);
void CMFEM_Ex10HyperelasticOperator_SetState(void *oper,
                                             const CMFEM_GridFunction *velocity,
                                             const CMFEM_GridFunction *position,
                                             CMFEM_Vector *state);
void CMFEM_Ex10HyperelasticOperator_SetFields(const void *oper,
                                              const CMFEM_Vector *state,
                                              CMFEM_GridFunction *velocity,
                                              CMFEM_GridFunction *position);
void CMFEM_Ex10HyperelasticOperator_Mult(void *oper,
                                         const CMFEM_Vector *state,
                                         CMFEM_Vector *dstate_dt);
void CMFEM_Ex10HyperelasticOperator_ImplicitSolve(void *oper,
                                                  double dt,
                                                  const CMFEM_Vector *state,
                                                  CMFEM_Vector *k);
double CMFEM_Ex10HyperelasticOperator_ElasticEnergy(void *oper,
                                                    const CMFEM_Vector *state);
double CMFEM_Ex10HyperelasticOperator_KineticEnergy(void *oper,
                                                    const CMFEM_Vector *state);
void CMFEM_Ex10HyperelasticOperator_GetElasticEnergyDensity(
   void *oper,
   const CMFEM_Vector *state,
   CMFEM_GridFunction *energy_density);

CMFEM_END_EXTERN_C

#endif
