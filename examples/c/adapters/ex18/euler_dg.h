#ifndef CMFEM_EULER_DG_H
#define CMFEM_EULER_DG_H

#include "../../../../c_api/common.h"

CMFEM_BEGIN_EXTERN_C

CMFEM_Mesh *CMFEM_Ex18_NewDefaultMesh(int problem);
void CMFEM_Ex18_ProjectInitialConditionGf(int problem,
                                          double specific_heat_ratio,
                                          double gas_constant,
                                          CMFEM_GridFunction *sol);
void *CMFEM_Ex18_NewDghcl(CMFEM_FiniteElementSpace *vfes,
                          int int_order_offset,
                          double specific_heat_ratio,
                          int preassemble_weak_div);
void CMFEM_Ex18_DeleteDghcl(void *context);
void CMFEM_Ex18_SetTimeDghcl(void *context, double time);
void CMFEM_Ex18_MultDghcl(void *context,
                          const CMFEM_Vector *x,
                          CMFEM_Vector *y);
double CMFEM_Ex18_GetMaxCharSpeedDghcl(void *context);
double CMFEM_Ex18_ComputeStableDtMeshDghclGf(const CMFEM_Mesh *mesh,
                                             int order,
                                             double cfl,
                                             void *context,
                                             const CMFEM_GridFunction *sol);
void CMFEM_Ex18_SaveStateMeshSfesGf(const CMFEM_Mesh *mesh,
                                    const CMFEM_FiniteElementSpace *scalar_fes,
                                    int num_equations,
                                    const CMFEM_GridFunction *sol,
                                    const char *mesh_path,
                                    const char *tag);
double CMFEM_Ex18_ComputeL2ErrorGf(int problem,
                                   double specific_heat_ratio,
                                   double gas_constant,
                                   const CMFEM_GridFunction *sol);
int CMFEM_Ex18_SendMomentumToGLVisMeshSfesDfesGf(
   const CMFEM_Mesh *mesh,
   const CMFEM_FiniteElementSpace *scalar_fes,
   const CMFEM_FiniteElementSpace *dfes,
   const CMFEM_GridFunction *sol,
   double time,
   int pause);

CMFEM_END_EXTERN_C

#endif
