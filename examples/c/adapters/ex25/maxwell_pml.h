#ifndef CMFEM_MAXWELL_PML_H
#define CMFEM_MAXWELL_PML_H

#include "../../../../c_api/common.h"

CMFEM_BEGIN_EXTERN_C

void *CMFEM_Ex25_NewPml(CMFEM_Mesh *mesh, int problem);
void CMFEM_Ex25_DeletePml(void *context);
void CMFEM_Ex25_SetAttributesPml(void *context, CMFEM_Mesh *mesh);
void CMFEM_Ex25_BuildEssentialBoundaryPml(void *context,
                                          const CMFEM_Mesh *mesh,
                                          CMFEM_ArrayInt *ess_bdr);
void CMFEM_Ex25_ProjectBoundaryPml(void *context,
                                   double mu,
                                   double epsilon,
                                   double omega,
                                   CMFEM_ComplexGridFunction *x,
                                   const CMFEM_ArrayInt *ess_bdr);
void CMFEM_Ex25_AssembleRhsPml(void *context,
                               double mu,
                               double epsilon,
                               double omega,
                               CMFEM_ComplexLinearForm *b);
void CMFEM_Ex25_AddOperatorIntegratorsPml(void *context,
                                          double mu,
                                          double epsilon,
                                          double omega,
                                          CMFEM_SesquilinearForm *a);
void CMFEM_Ex25_AddPreconditionerIntegratorsPml(void *context,
                                                double mu,
                                                double epsilon,
                                                double omega,
                                                CMFEM_BilinearForm *prec);
void CMFEM_Ex25_ComputeErrorsPml(void *context,
                                 double mu,
                                 double epsilon,
                                 double omega,
                                 const CMFEM_ComplexGridFunction *x,
                                 int order,
                                 double *rel_re,
                                 double *rel_im,
                                 double *total);

CMFEM_END_EXTERN_C

#endif
