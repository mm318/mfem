#ifndef CMFEM_MESH_H
#define CMFEM_MESH_H

#include "common.h"

CMFEM_STORAGE(Mesh, 5712);

CMFEM_BEGIN_EXTERN_C

CMFEM_Mesh CMFEM_Mesh_Construct(void);
CMFEM_Mesh *CMFEM_Mesh_New(void);
CMFEM_Mesh *CMFEM_Mesh_NewFile(const char *mesh_file, int generate_edges,
                               int refine);
CMFEM_Mesh CMFEM_Mesh_Copy(const CMFEM_Mesh *mesh);
CMFEM_Mesh *CMFEM_Mesh_NewCopy(const CMFEM_Mesh *mesh);
void CMFEM_Mesh_Destroy(CMFEM_Mesh *mesh);
void CMFEM_Mesh_Delete(CMFEM_Mesh *mesh);
int CMFEM_Mesh_Dimension(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_SpaceDimension(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_GetNE(const CMFEM_Mesh *mesh);
void CMFEM_Mesh_UniformRefinement(CMFEM_Mesh *mesh);
int CMFEM_Mesh_BoundaryAttributesSize(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_BoundaryAttributesMax(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_AttributesMax(const CMFEM_Mesh *mesh);
void CMFEM_Mesh_MarkExternalBoundaries(const CMFEM_Mesh *mesh,
                                       CMFEM_ArrayInt *ess_bdr);
int CMFEM_Mesh_HasNURBSext(const CMFEM_Mesh *mesh);
void CMFEM_Mesh_DegreeElevate(CMFEM_Mesh *mesh, int order, int max_degree);
void CMFEM_Mesh_EnsureNCMesh(CMFEM_Mesh *mesh, int simplices_nonconforming);
int CMFEM_Mesh_HasNodes(const CMFEM_Mesh *mesh);
CMFEM_FiniteElementSpace *CMFEM_Mesh_GetNodesFESpace(CMFEM_Mesh *mesh);
const char *CMFEM_Mesh_GetNodesOwnFECName(CMFEM_Mesh *mesh);
void CMFEM_Mesh_SetNodalFESpace(CMFEM_Mesh *mesh,
                                CMFEM_FiniteElementSpace *fespace);
void CMFEM_Mesh_AddDisplacementToNodes(CMFEM_Mesh *mesh,
                                       const CMFEM_GridFunction *displacement);
void CMFEM_Mesh_Save(const CMFEM_Mesh *mesh, const char *path, int precision);
void CMFEM_Mesh_Print(const CMFEM_Mesh *mesh, const char *path, int precision);

CMFEM_END_EXTERN_C

#endif
