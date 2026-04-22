#include "common.hpp"

extern "C" {

   CMFEM_SparseMatrix *CMFEM_BuildNormalConstraints(
      CMFEM_FiniteElementSpace *fespace,
      CMFEM_ArrayInt *constrained_att,
      CMFEM_ArrayInt *constraint_rowstarts)
   {
      return reinterpret_cast<CMFEM_SparseMatrix *>(
                mfem::BuildNormalConstraints(
                   *cmfem::As<mfem::FiniteElementSpace>(fespace),
                   cmfem::ArrayIntRef(constrained_att),
                   cmfem::ArrayIntRef(constraint_rowstarts)));
   }

} // extern "C"
