#include "common.hpp"

#include <fstream>

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_Mesh, mfem::Mesh);

mfem::AttributeSets &AttributeSetsRef(mfem::Mesh *mesh, int boundary)
{
   return boundary != 0 ? mesh->bdr_attribute_sets : mesh->attribute_sets;
}

const mfem::AttributeSets &AttributeSetsRef(const mfem::Mesh *mesh,
                                            int boundary)
{
   return boundary != 0 ? mesh->bdr_attribute_sets : mesh->attribute_sets;
}

void SnapNodesToUnitSphere(mfem::Mesh &mesh)
{
   mfem::GridFunction &nodes = *mesh.GetNodes();
   mfem::Vector node(mesh.SpaceDimension());
   for (int i = 0; i < nodes.FESpace()->GetNDofs(); i++)
   {
      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
         node(d) = nodes(nodes.FESpace()->DofToVDof(i, d));
      }

      node /= node.Norml2();

      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
         nodes(nodes.FESpace()->DofToVDof(i, d)) = node(d);
      }
   }
   if (mesh.Nonconforming())
   {
      mfem::Vector tnodes;
      nodes.GetTrueDofs(tnodes);
      nodes.SetFromTrueDofs(tnodes);
   }
}

} // namespace

extern "C" {

   CMFEM_Mesh CMFEM_Mesh_Construct(void)
   {
      alignas(mfem::Mesh) CMFEM_Mesh mesh;
      new (cmfem::As<mfem::Mesh>(&mesh)) mfem::Mesh();
      return mesh;
   }

   CMFEM_Mesh *CMFEM_Mesh_New(void)
   {
      return reinterpret_cast<CMFEM_Mesh *>(new mfem::Mesh());
   }

   CMFEM_Mesh *CMFEM_Mesh_NewDimensionVerticesElementsBoundarySpace(
      int dim,
      int num_vertices,
      int num_elements,
      int num_bdr_elements,
      int space_dim)
   {
      return reinterpret_cast<CMFEM_Mesh *>(
                new mfem::Mesh(dim, num_vertices, num_elements,
                               num_bdr_elements, space_dim));
   }

   CMFEM_Mesh *CMFEM_Mesh_NewFile(const char *mesh_file, int generate_edges,
                                  int refine)
   {
      return reinterpret_cast<CMFEM_Mesh *>(
                new mfem::Mesh(mesh_file, generate_edges, refine));
   }

   CMFEM_Mesh CMFEM_Mesh_Copy(const CMFEM_Mesh *mesh)
   {
      alignas(mfem::Mesh) CMFEM_Mesh copy;
      new (cmfem::As<mfem::Mesh>(&copy)) mfem::Mesh(*cmfem::As<const mfem::Mesh>
                                                    (mesh));
      return copy;
   }

   CMFEM_Mesh *CMFEM_Mesh_NewCopy(const CMFEM_Mesh *mesh)
   {
      return reinterpret_cast<CMFEM_Mesh *>(
                new mfem::Mesh(*cmfem::As<const mfem::Mesh>(mesh)));
   }

   void CMFEM_Mesh_Destroy(CMFEM_Mesh *mesh)
   {
      cmfem::As<mfem::Mesh>(mesh)->~Mesh();
   }

   void CMFEM_Mesh_Delete(CMFEM_Mesh *mesh)
   {
      delete cmfem::As<mfem::Mesh>(mesh);
   }

   int CMFEM_Mesh_Dimension(const CMFEM_Mesh *mesh)
   {
      return cmfem::As<const mfem::Mesh>(mesh)->Dimension();
   }

   int CMFEM_Mesh_SpaceDimension(const CMFEM_Mesh *mesh)
   {
      return cmfem::As<const mfem::Mesh>(mesh)->SpaceDimension();
   }

   int CMFEM_Mesh_GetNE(const CMFEM_Mesh *mesh)
   {
      return cmfem::As<const mfem::Mesh>(mesh)->GetNE();
   }

   int CMFEM_Mesh_GetElementType(const CMFEM_Mesh *mesh, int index)
   {
      return static_cast<int>(cmfem::As<const mfem::Mesh>(mesh)->GetElementType(
                                 index));
   }

   void CMFEM_Mesh_AddVertex(CMFEM_Mesh *mesh, const double *coords)
   {
      cmfem::As<mfem::Mesh>(mesh)->AddVertex(coords);
   }

   void CMFEM_Mesh_AddTriangle(CMFEM_Mesh *mesh, const int *indices,
                               int attribute)
   {
      cmfem::As<mfem::Mesh>(mesh)->AddTriangle(indices, attribute);
   }

   void CMFEM_Mesh_AddQuad(CMFEM_Mesh *mesh, const int *indices, int attribute)
   {
      cmfem::As<mfem::Mesh>(mesh)->AddQuad(indices, attribute);
   }

   void CMFEM_Mesh_FinalizeTriMesh(CMFEM_Mesh *mesh, int generate_edges,
                                   int refine, int fix_orientation)
   {
      cmfem::As<mfem::Mesh>(mesh)->FinalizeTriMesh(generate_edges, refine,
                                                   fix_orientation != 0);
   }

   void CMFEM_Mesh_FinalizeQuadMesh(CMFEM_Mesh *mesh, int generate_edges,
                                    int refine, int fix_orientation)
   {
      cmfem::As<mfem::Mesh>(mesh)->FinalizeQuadMesh(generate_edges, refine,
                                                    fix_orientation != 0);
   }

   void CMFEM_Mesh_UniformRefinement(CMFEM_Mesh *mesh)
   {
      cmfem::As<mfem::Mesh>(mesh)->UniformRefinement();
   }

   void CMFEM_Mesh_Finalize(CMFEM_Mesh *mesh, int refine, int fix_orientation)
   {
      cmfem::As<mfem::Mesh>(mesh)->Finalize(refine != 0, fix_orientation != 0);
   }

   void CMFEM_Mesh_RefineAtVertex3(CMFEM_Mesh *mesh, double x, double y,
                                   double z)
   {
      cmfem::As<mfem::Mesh>(mesh)->RefineAtVertex(mfem::Vertex(
                                                     static_cast<mfem::real_t>(x),
                                                     static_cast<mfem::real_t>(y),
                                                     static_cast<mfem::real_t>(z)));
   }

   void CMFEM_Mesh_RandomRefinement(CMFEM_Mesh *mesh, double probability)
   {
      cmfem::As<mfem::Mesh>(mesh)->RandomRefinement(
         static_cast<mfem::real_t>(probability));
   }

   int CMFEM_Mesh_BoundaryAttributesSize(const CMFEM_Mesh *mesh)
   {
      return cmfem::As<const mfem::Mesh>(mesh)->bdr_attributes.Size();
   }

   int CMFEM_Mesh_BoundaryAttributesMax(const CMFEM_Mesh *mesh)
   {
      return cmfem::As<const mfem::Mesh>(mesh)->bdr_attributes.Max();
   }

   int CMFEM_Mesh_AttributesMax(const CMFEM_Mesh *mesh)
   {
      return cmfem::As<const mfem::Mesh>(mesh)->attributes.Max();
   }

   void CMFEM_Mesh_MarkExternalBoundaries(const CMFEM_Mesh *mesh,
                                          CMFEM_ArrayInt *ess_bdr)
   {
      cmfem::As<const mfem::Mesh>(mesh)->MarkExternalBoundaries(cmfem::ArrayIntRef(
                                                                   ess_bdr));
   }

   int CMFEM_Mesh_HasNURBSext(const CMFEM_Mesh *mesh)
   {
      return cmfem::As<const mfem::Mesh>(mesh)->NURBSext != nullptr;
   }

   void CMFEM_Mesh_DegreeElevate(CMFEM_Mesh *mesh, int order, int max_degree)
   {
      cmfem::As<mfem::Mesh>(mesh)->DegreeElevate(order, max_degree);
   }

   void CMFEM_Mesh_EnsureNCMesh(CMFEM_Mesh *mesh, int simplices_nonconforming)
   {
      cmfem::As<mfem::Mesh>(mesh)->EnsureNCMesh(simplices_nonconforming != 0);
   }

   void CMFEM_Mesh_SetCurvature(CMFEM_Mesh *mesh, int order)
   {
      cmfem::As<mfem::Mesh>(mesh)->SetCurvature(order);
   }

   int CMFEM_Mesh_HasNodes(const CMFEM_Mesh *mesh)
   {
      return cmfem::As<const mfem::Mesh>(mesh)->GetNodes() != nullptr;
   }

   CMFEM_FiniteElementSpace *CMFEM_Mesh_GetNodesFESpace(CMFEM_Mesh *mesh)
   {
      mfem::GridFunction *nodes = cmfem::As<mfem::Mesh>(mesh)->GetNodes();
      return nodes ? reinterpret_cast<CMFEM_FiniteElementSpace *>
             (nodes->FESpace()) : nullptr;
   }

   const char *CMFEM_Mesh_GetNodesOwnFECName(CMFEM_Mesh *mesh)
   {
      mfem::GridFunction *nodes = cmfem::As<mfem::Mesh>(mesh)->GetNodes();
      return (nodes && nodes->OwnFEC()) ? nodes->OwnFEC()->Name() : "";
   }

   void CMFEM_Mesh_SetNodalFESpace(CMFEM_Mesh *mesh,
                                   CMFEM_FiniteElementSpace *fespace)
   {
      cmfem::As<mfem::Mesh>(mesh)->SetNodalFESpace(
         cmfem::As<mfem::FiniteElementSpace>(fespace));
   }

   void CMFEM_Mesh_SnapNodesToUnitSphere(CMFEM_Mesh *mesh)
   {
      SnapNodesToUnitSphere(*cmfem::As<mfem::Mesh>(mesh));
   }

   void CMFEM_Mesh_AddDisplacementToNodes(CMFEM_Mesh *mesh,
                                          const CMFEM_GridFunction *displacement)
   {
      mfem::GridFunction *nodes = cmfem::As<mfem::Mesh>(mesh)->GetNodes();
      if (nodes)
      {
         *nodes += *cmfem::As<const mfem::GridFunction>(displacement);
      }
   }

   void CMFEM_Mesh_PrintAttributeSetNames(const CMFEM_Mesh *mesh, int boundary,
                                          const char *label)
   {
      const auto &sets = AttributeSetsRef(cmfem::As<const mfem::Mesh>(mesh),
                                          boundary);
      const auto names = sets.GetAttributeSetNames();
      if (label && label[0] != '\0')
      {
         mfem::out << label << ": ";
      }
      for (const auto &name : names)
      {
         mfem::out << " \"" << name << "\"";
      }
      mfem::out << '\n';
   }

   int CMFEM_Mesh_AttributeSetExists(const CMFEM_Mesh *mesh, int boundary,
                                     const char *name)
   {
      return AttributeSetsRef(cmfem::As<const mfem::Mesh>(mesh), boundary)
             .AttributeSetExists(name) ? 1 : 0;
   }

   void CMFEM_Mesh_SetAttributeSetFromName(CMFEM_Mesh *mesh, int boundary,
                                           const char *target_name,
                                           const char *source_name)
   {
      auto &sets = AttributeSetsRef(cmfem::As<mfem::Mesh>(mesh), boundary);
      sets.SetAttributeSet(target_name, sets.GetAttributeSet(source_name));
   }

   void CMFEM_Mesh_AddToAttributeSetFromName(CMFEM_Mesh *mesh, int boundary,
                                             const char *target_name,
                                             const char *source_name)
   {
      auto &sets = AttributeSetsRef(cmfem::As<mfem::Mesh>(mesh), boundary);
      sets.AddToAttributeSet(target_name, sets.GetAttributeSet(source_name));
   }

   void CMFEM_Mesh_GetAttributeSetMarker(const CMFEM_Mesh *mesh, int boundary,
                                         const char *name,
                                         CMFEM_ArrayInt *marker)
   {
      const auto &sets = AttributeSetsRef(cmfem::As<const mfem::Mesh>(mesh),
                                          boundary);
      cmfem::ArrayIntRef(marker) = sets.GetAttributeSetMarker(name);
   }

   void CMFEM_Mesh_Save(const CMFEM_Mesh *mesh, const char *path, int precision)
   {
      cmfem::As<const mfem::Mesh>(mesh)->Save(path, precision);
   }

   void CMFEM_Mesh_Print(const CMFEM_Mesh *mesh, const char *path, int precision)
   {
      std::ofstream out(path);
      out.precision(precision);
      cmfem::As<const mfem::Mesh>(mesh)->Print(out);
   }

} // extern "C"
