#ifndef CMFEM_FINITE_ELEMENT_SPACE_H
#define CMFEM_FINITE_ELEMENT_SPACE_H

#include "common.h"

CMFEM_STORAGE(FiniteElementSpace, 1400);

CMFEM_BEGIN_EXTERN_C

CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshH1(CMFEM_Mesh *mesh,
                                                             CMFEM_H1FeCollection *fec);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshH1VDim(
   CMFEM_Mesh *mesh, CMFEM_H1FeCollection *fec, int vdim);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshNd(CMFEM_Mesh *mesh,
                                                             CMFEM_NdFeCollection *fec);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshRt(CMFEM_Mesh *mesh,
                                                             CMFEM_RtFeCollection *fec);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshL2Vdim(
   CMFEM_Mesh *mesh,
   CMFEM_L2FeCollection *fec,
   int vdim);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshDg(
   CMFEM_Mesh *mesh, CMFEM_DgFeCollection *fec);
CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshDgInterface(
   CMFEM_Mesh *mesh, CMFEM_DgInterfaceFeCollection *fec);
void CMFEM_FiniteElementSpace_Delete(CMFEM_FiniteElementSpace *fespace);
int CMFEM_FiniteElementSpace_GetVSize(const CMFEM_FiniteElementSpace *fespace);
int CMFEM_FiniteElementSpace_GetTrueVSize(const CMFEM_FiniteElementSpace
                                          *fespace);
void CMFEM_FiniteElementSpace_GetBoundaryTrueDofs(const CMFEM_FiniteElementSpace
                                                  *fespace, CMFEM_ArrayInt *boundary_dofs);
void CMFEM_FiniteElementSpace_GetEssentialTrueDofs(const
                                                   CMFEM_FiniteElementSpace *fespace, const CMFEM_ArrayInt *ess_bdr,
                                                   CMFEM_ArrayInt *ess_tdof_list);
void CMFEM_FiniteElementSpace_Update(CMFEM_FiniteElementSpace *fespace);
void CMFEM_FiniteElementSpace_UpdatesFinished(CMFEM_FiniteElementSpace
                                              *fespace);

CMFEM_END_EXTERN_C

#endif
