#ifndef CMFEM_MESH_H
#define CMFEM_MESH_H

#include "common.h"

CMFEM_STORAGE(Mesh, 5712);

typedef void (*CMFEM_MeshTransformCallback)(const CMFEM_Vector *x,
                                            CMFEM_Vector *value,
                                            void *context);

enum CMFEM_ElementType
{
   CMFEM_ELEMENT_POINT = 0,
   CMFEM_ELEMENT_SEGMENT = 1,
   CMFEM_ELEMENT_TRIANGLE = 2,
   CMFEM_ELEMENT_QUADRILATERAL = 3,
   CMFEM_ELEMENT_TETRAHEDRON = 4,
   CMFEM_ELEMENT_HEXAHEDRON = 5,
   CMFEM_ELEMENT_WEDGE = 6,
   CMFEM_ELEMENT_PYRAMID = 7,
};

CMFEM_BEGIN_EXTERN_C

CMFEM_Mesh CMFEM_Mesh_Construct(void);
CMFEM_Mesh *CMFEM_Mesh_New(void);
CMFEM_Mesh *CMFEM_Mesh_NewDimensionVerticesElementsBoundarySpace(int dim,
                                                                 int num_vertices,
                                                                 int num_elements,
                                                                 int num_bdr_elements,
                                                                 int space_dim);
CMFEM_Mesh *CMFEM_Mesh_NewFile(const char *mesh_file, int generate_edges,
                               int refine);
CMFEM_Mesh CMFEM_Mesh_Copy(const CMFEM_Mesh *mesh);
CMFEM_Mesh *CMFEM_Mesh_NewCopy(const CMFEM_Mesh *mesh);
void CMFEM_Mesh_Destroy(CMFEM_Mesh *mesh);
void CMFEM_Mesh_Delete(CMFEM_Mesh *mesh);
int CMFEM_Mesh_Dimension(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_SpaceDimension(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_MeshGenerator(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_GetNE(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_GetElementType(const CMFEM_Mesh *mesh, int index);
void CMFEM_Mesh_AddVertex(CMFEM_Mesh *mesh, const double *coords);
void CMFEM_Mesh_AddTriangle(CMFEM_Mesh *mesh, const int *indices,
                            int attribute);
void CMFEM_Mesh_AddQuad(CMFEM_Mesh *mesh, const int *indices, int attribute);
void CMFEM_Mesh_AddBdrSegment(CMFEM_Mesh *mesh, const int *indices,
                              int attribute);
void CMFEM_Mesh_FinalizeTopology(CMFEM_Mesh *mesh);
void CMFEM_Mesh_FinalizeTriMesh(CMFEM_Mesh *mesh, int generate_edges,
                                int refine, int fix_orientation);
void CMFEM_Mesh_FinalizeQuadMesh(CMFEM_Mesh *mesh, int generate_edges,
                                 int refine, int fix_orientation);
void CMFEM_Mesh_UniformRefinement(CMFEM_Mesh *mesh);
void CMFEM_Mesh_Finalize(CMFEM_Mesh *mesh, int refine, int fix_orientation);
void CMFEM_Mesh_RefineAtVertex3(CMFEM_Mesh *mesh, double x, double y, double z);
void CMFEM_Mesh_RandomRefinement(CMFEM_Mesh *mesh, double probability);
int CMFEM_Mesh_BoundaryAttributesSize(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_BoundaryAttributesMax(const CMFEM_Mesh *mesh);
int CMFEM_Mesh_AttributesMax(const CMFEM_Mesh *mesh);
void CMFEM_Mesh_MarkExternalBoundaries(const CMFEM_Mesh *mesh,
                                       CMFEM_ArrayInt *ess_bdr);
int CMFEM_Mesh_HasNURBSext(const CMFEM_Mesh *mesh);
void CMFEM_Mesh_DegreeElevate(CMFEM_Mesh *mesh, int order, int max_degree);
void CMFEM_Mesh_EnsureNCMesh(CMFEM_Mesh *mesh, int simplices_nonconforming);
void CMFEM_Mesh_SetCurvature(CMFEM_Mesh *mesh, int order);
int CMFEM_Mesh_HasNodes(const CMFEM_Mesh *mesh);
CMFEM_FiniteElementSpace *CMFEM_Mesh_GetNodesFESpace(CMFEM_Mesh *mesh);
const char *CMFEM_Mesh_GetNodesOwnFECName(CMFEM_Mesh *mesh);
void CMFEM_Mesh_GetNodesGf(const CMFEM_Mesh *mesh,
                           CMFEM_GridFunction *nodes);
void CMFEM_Mesh_GetBoundingBox(const CMFEM_Mesh *mesh,
                               CMFEM_Vector *min,
                               CMFEM_Vector *max,
                               int ref);
void CMFEM_Mesh_SetNodalFESpace(CMFEM_Mesh *mesh,
                                CMFEM_FiniteElementSpace *fespace);
void CMFEM_Mesh_Transform(CMFEM_Mesh *mesh,
                          CMFEM_MeshTransformCallback callback,
                          void *context);
void CMFEM_Mesh_SnapNodesToUnitSphere(CMFEM_Mesh *mesh);
void CMFEM_Mesh_AddDisplacementToNodes(CMFEM_Mesh *mesh,
                                       const CMFEM_GridFunction *displacement);
void CMFEM_Mesh_PrintAttributeSetNames(const CMFEM_Mesh *mesh, int boundary,
                                       const char *label);
int CMFEM_Mesh_AttributeSetExists(const CMFEM_Mesh *mesh, int boundary,
                                  const char *name);
void CMFEM_Mesh_SetAttributeSetFromName(CMFEM_Mesh *mesh, int boundary,
                                        const char *target_name,
                                        const char *source_name);
void CMFEM_Mesh_AddToAttributeSetFromName(CMFEM_Mesh *mesh, int boundary,
                                          const char *target_name,
                                          const char *source_name);
void CMFEM_Mesh_GetAttributeSetMarker(const CMFEM_Mesh *mesh, int boundary,
                                      const char *name,
                                      CMFEM_ArrayInt *marker);
void CMFEM_Mesh_Save(const CMFEM_Mesh *mesh, const char *path, int precision);
void CMFEM_Mesh_Print(const CMFEM_Mesh *mesh, const char *path, int precision);

CMFEM_END_EXTERN_C

#endif
