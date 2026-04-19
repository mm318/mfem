#ifndef CMFEM_MAGNETOSTATICS_SUBMESH_H
#define CMFEM_MAGNETOSTATICS_SUBMESH_H

#include "common.h"

CMFEM_BEGIN_EXTERN_C

void CMFEM_Ex34_PrepareCurrentDensityRt(CMFEM_Mesh *mesh,
                                        int order,
                                        int pa,
                                        int visualization,
                                        const CMFEM_ArrayInt *submesh_elems,
                                        const CMFEM_ArrayInt *phi0_attr,
                                        const CMFEM_ArrayInt *phi1_attr,
                                        const CMFEM_ArrayInt *jn_zero_attr,
                                        CMFEM_GridFunction *j_full);
void CMFEM_Ex34_AssembleCurrentDensityLf(CMFEM_FiniteElementSpace *fespace_nd,
                                         const CMFEM_GridFunction *j_full,
                                         CMFEM_LinearForm *linear_form);

CMFEM_END_EXTERN_C

#endif
