#ifndef CMFEM_FINITE_ELEMENT_SPACE_H
#define CMFEM_FINITE_ELEMENT_SPACE_H

#include "common.h"

CMFEM_STORAGE(FiniteElementSpace, 1400);

CMFEM_BEGIN_EXTERN_C

CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshH1(CMFEM_Mesh *mesh,
                                                             CMFEM_H1_FECollection *fec);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshH1VDim(
   CMFEM_Mesh *mesh, CMFEM_H1_FECollection *fec, int vdim);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshND(CMFEM_Mesh *mesh,
                                                             CMFEM_ND_FECollection *fec);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshRT(CMFEM_Mesh *mesh,
                                                             CMFEM_RT_FECollection *fec);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshDG(
   CMFEM_Mesh *mesh, CMFEM_DG_FECollection *fec);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshDGInterface(
   CMFEM_Mesh *mesh, CMFEM_DG_Interface_FECollection *fec);
void CMFEM_FiniteElementSpace_Delete(CMFEM_FiniteElementSpace *fespace);
int CMFEM_FiniteElementSpace_GetTrueVSize(const CMFEM_FiniteElementSpace
                                          *fespace);
void CMFEM_FiniteElementSpace_GetBoundaryTrueDofs(const CMFEM_FiniteElementSpace
                                                  *fespace, CMFEM_ArrayInt *boundary_dofs);
void CMFEM_FiniteElementSpace_GetEssentialTrueDofs(const
                                                   CMFEM_FiniteElementSpace *fespace, const CMFEM_ArrayInt *ess_bdr,
                                                   CMFEM_ArrayInt *ess_tdof_list);

CMFEM_END_EXTERN_C

#endif
