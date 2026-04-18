#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_FiniteElementSpace, mfem::FiniteElementSpace);

} // namespace

extern "C" {

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshH1(CMFEM_Mesh *mesh,
                                                                CMFEM_H1FeCollection *fec)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::H1_FECollection>(fec)));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshH1VDim(
      CMFEM_Mesh *mesh,
      CMFEM_H1FeCollection *fec,
      int vdim)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::H1_FECollection>(fec),
                                             vdim));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshNd(CMFEM_Mesh *mesh,
                                                                CMFEM_NdFeCollection *fec)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::ND_FECollection>(fec)));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshNdR1d(
      CMFEM_Mesh *mesh,
      CMFEM_NdR1dFeCollection *fec)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::ND_R1D_FECollection>(fec)));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshNdR2d(
      CMFEM_Mesh *mesh,
      CMFEM_NdR2dFeCollection *fec)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::ND_R2D_FECollection>(fec)));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshRt(CMFEM_Mesh *mesh,
                                                                CMFEM_RtFeCollection *fec)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::RT_FECollection>(fec)));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshRtTrace(
      CMFEM_Mesh *mesh,
      CMFEM_RtTraceFeCollection *fec)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(
                   cmfem::As<mfem::Mesh>(mesh),
                   cmfem::As<mfem::RT_Trace_FECollection>(fec)));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshL2Vdim(
      CMFEM_Mesh *mesh,
      CMFEM_L2FeCollection *fec,
      int vdim)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::L2_FECollection>(fec),
                                             vdim));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshDg(CMFEM_Mesh *mesh,
                                                                CMFEM_DgFeCollection *fec)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::DG_FECollection>(fec)));
   }

   CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshDgVDim(
      CMFEM_Mesh *mesh,
      CMFEM_DgFeCollection *fec,
      int vdim)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::DG_FECollection>(fec),
                                             vdim));
   }

   CMFEM_FiniteElementSpace *
   CMFEM_FiniteElementSpace_NewMeshDgInterface(CMFEM_Mesh *mesh,
                                               CMFEM_DgInterfaceFeCollection *fec)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                             cmfem::As<mfem::DG_Interface_FECollection>(fec)));
   }

   void CMFEM_FiniteElementSpace_Delete(CMFEM_FiniteElementSpace *fespace)
   {
      delete cmfem::As<mfem::FiniteElementSpace>(fespace);
   }

   int CMFEM_FiniteElementSpace_GetVSize(const CMFEM_FiniteElementSpace
                                         *fespace)
   {
      return cmfem::As<const mfem::FiniteElementSpace>(fespace)->GetVSize();
   }

   int CMFEM_FiniteElementSpace_GetTrueVSize(const CMFEM_FiniteElementSpace
                                             *fespace)
   {
      return cmfem::As<const mfem::FiniteElementSpace>(fespace)->GetTrueVSize();
   }

   const CMFEM_SparseMatrix *CMFEM_FiniteElementSpace_GetRestrictionMatrixSm(
      const CMFEM_FiniteElementSpace *fespace)
   {
      return reinterpret_cast<const CMFEM_SparseMatrix *>(
                cmfem::As<const mfem::FiniteElementSpace>(fespace)
                ->GetRestrictionMatrix());
   }

   void CMFEM_FiniteElementSpace_GetBoundaryTrueDofs(const CMFEM_FiniteElementSpace
                                                     *fespace,
                                                     CMFEM_ArrayInt *boundary_dofs)
   {
      const_cast<mfem::FiniteElementSpace *>(
         cmfem::As<const mfem::FiniteElementSpace>(fespace))->GetBoundaryTrueDofs(
            cmfem::ArrayIntRef(boundary_dofs));
   }

   void CMFEM_FiniteElementSpace_GetEssentialTrueDofs(const
                                                      CMFEM_FiniteElementSpace *fespace,
                                                      const CMFEM_ArrayInt *ess_bdr,
                                                      CMFEM_ArrayInt *ess_tdof_list)
   {
      cmfem::As<const mfem::FiniteElementSpace>(fespace)->GetEssentialTrueDofs(
         cmfem::ArrayIntRef(ess_bdr),
         cmfem::ArrayIntRef(ess_tdof_list));
   }

   void CMFEM_FiniteElementSpace_Update(CMFEM_FiniteElementSpace *fespace)
   {
      cmfem::As<mfem::FiniteElementSpace>(fespace)->Update();
   }

   void CMFEM_FiniteElementSpace_UpdatesFinished(CMFEM_FiniteElementSpace *fespace)
   {
      cmfem::As<mfem::FiniteElementSpace>(fespace)->UpdatesFinished();
   }

} // extern "C"
