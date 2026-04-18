#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_CGSolver, mfem::CGSolver);

} // namespace

extern "C" {

   CMFEM_CGSolver *CMFEM_CGSolver_New(void)
   {
      return reinterpret_cast<CMFEM_CGSolver *>(new mfem::CGSolver());
   }

   void CMFEM_CGSolver_Delete(CMFEM_CGSolver *solver)
   {
      delete cmfem::As<mfem::CGSolver>(solver);
   }

   void CMFEM_CGSolver_SetOperatorSm(CMFEM_CGSolver *solver,
                                     const CMFEM_SparseMatrix *matrix)
   {
      cmfem::As<mfem::CGSolver>(solver)->SetOperator(
         *cmfem::As<const mfem::SparseMatrix>(matrix));
   }

   void CMFEM_CGSolver_SetPrintLevel(CMFEM_CGSolver *solver, int print_level)
   {
      cmfem::As<mfem::CGSolver>(solver)->SetPrintLevel(print_level);
   }

   void CMFEM_CGSolver_SetRelTol(CMFEM_CGSolver *solver, double rel_tol)
   {
      cmfem::As<mfem::CGSolver>(solver)->SetRelTol(
         static_cast<mfem::real_t>(rel_tol));
   }

   void CMFEM_CGSolver_SetMaxIter(CMFEM_CGSolver *solver, int max_iter)
   {
      cmfem::As<mfem::CGSolver>(solver)->SetMaxIter(max_iter);
   }

   void CMFEM_CGSolver_SetIterativeMode(CMFEM_CGSolver *solver,
                                        int iterative_mode)
   {
      cmfem::As<mfem::CGSolver>(solver)->iterative_mode = iterative_mode != 0;
   }

} // extern "C"
