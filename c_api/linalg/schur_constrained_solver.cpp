#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_SchurConstrainedSolver, mfem::SchurConstrainedSolver);

} // namespace

extern "C" {

   CMFEM_SchurConstrainedSolver *CMFEM_SchurConstrainedSolver_NewSmSmGs(
      CMFEM_SparseMatrix *primal_operator,
      CMFEM_SparseMatrix *constraint_operator,
      CMFEM_GSSmoother *primal_pc)
   {
      return reinterpret_cast<CMFEM_SchurConstrainedSolver *>(
                new mfem::SchurConstrainedSolver(
                   *cmfem::As<mfem::SparseMatrix>(primal_operator),
                   *cmfem::As<mfem::SparseMatrix>(constraint_operator),
                   *cmfem::As<mfem::GSSmoother>(primal_pc)));
   }

   void CMFEM_SchurConstrainedSolver_Delete(
      CMFEM_SchurConstrainedSolver *solver)
   {
      delete cmfem::As<mfem::SchurConstrainedSolver>(solver);
   }

   void CMFEM_SchurConstrainedSolver_SetRelTol(
      CMFEM_SchurConstrainedSolver *solver,
      double rel_tol)
   {
      cmfem::As<mfem::SchurConstrainedSolver>(solver)->SetRelTol(
         static_cast<mfem::real_t>(rel_tol));
   }

   void CMFEM_SchurConstrainedSolver_SetMaxIter(
      CMFEM_SchurConstrainedSolver *solver,
      int max_iter)
   {
      cmfem::As<mfem::SchurConstrainedSolver>(solver)->SetMaxIter(max_iter);
   }

   void CMFEM_SchurConstrainedSolver_SetPrintLevel(
      CMFEM_SchurConstrainedSolver *solver,
      int print_level)
   {
      cmfem::As<mfem::SchurConstrainedSolver>(solver)->SetPrintLevel(
         print_level);
   }

   void CMFEM_SchurConstrainedSolver_Mult(
      const CMFEM_SchurConstrainedSolver *solver,
      const CMFEM_Vector *rhs,
      CMFEM_Vector *solution)
   {
      cmfem::As<const mfem::SchurConstrainedSolver>(solver)->Mult(
         cmfem::VectorRef(rhs),
         cmfem::VectorRef(solution));
   }

} // extern "C"
