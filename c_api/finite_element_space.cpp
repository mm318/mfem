#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_FiniteElementSpace, mfem::FiniteElementSpace);

} // namespace

extern "C" {

CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshH1(CMFEM_Mesh *mesh,
                                                             CMFEM_H1_FECollection *fec)
{
   return reinterpret_cast<CMFEM_FiniteElementSpace *>(
      new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                   cmfem::As<mfem::H1_FECollection>(fec)));
}

CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshH1VDim(CMFEM_Mesh *mesh,
                                                                 CMFEM_H1_FECollection *fec,
                                                                 int vdim)
{
   return reinterpret_cast<CMFEM_FiniteElementSpace *>(
      new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                   cmfem::As<mfem::H1_FECollection>(fec),
                                   vdim));
}

CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshND(CMFEM_Mesh *mesh,
                                                             CMFEM_ND_FECollection *fec)
{
   return reinterpret_cast<CMFEM_FiniteElementSpace *>(
      new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                   cmfem::As<mfem::ND_FECollection>(fec)));
}

CMFEM_FiniteElementSpace *CMFEM_FiniteElementSpace_NewMeshRT(CMFEM_Mesh *mesh,
                                                             CMFEM_RT_FECollection *fec)
{
   return reinterpret_cast<CMFEM_FiniteElementSpace *>(
      new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                   cmfem::As<mfem::RT_FECollection>(fec)));
}

CMFEM_FiniteElementSpace *
CMFEM_FiniteElementSpace_NewMeshDGInterface(CMFEM_Mesh *mesh,
                                            CMFEM_DG_Interface_FECollection *fec)
{
   return reinterpret_cast<CMFEM_FiniteElementSpace *>(
      new mfem::FiniteElementSpace(cmfem::As<mfem::Mesh>(mesh),
                                   cmfem::As<mfem::DG_Interface_FECollection>(fec)));
}

void CMFEM_FiniteElementSpace_Delete(CMFEM_FiniteElementSpace *fespace)
{
   delete cmfem::As<mfem::FiniteElementSpace>(fespace);
}

int CMFEM_FiniteElementSpace_GetTrueVSize(const CMFEM_FiniteElementSpace *fespace)
{
   return cmfem::As<const mfem::FiniteElementSpace>(fespace)->GetTrueVSize();
}

void CMFEM_FiniteElementSpace_GetBoundaryTrueDofs(const CMFEM_FiniteElementSpace *fespace,
                                                  CMFEM_ArrayInt *boundary_dofs)
{
   const_cast<mfem::FiniteElementSpace *>(
      cmfem::As<const mfem::FiniteElementSpace>(fespace))->GetBoundaryTrueDofs(
      cmfem::ArrayIntRef(boundary_dofs));
}

void CMFEM_FiniteElementSpace_GetEssentialTrueDofs(const CMFEM_FiniteElementSpace *fespace,
                                                   const CMFEM_ArrayInt *ess_bdr,
                                                   CMFEM_ArrayInt *ess_tdof_list)
{
   cmfem::As<const mfem::FiniteElementSpace>(fespace)->GetEssentialTrueDofs(
      cmfem::ArrayIntRef(ess_bdr),
      cmfem::ArrayIntRef(ess_tdof_list));
}

} // extern "C"
