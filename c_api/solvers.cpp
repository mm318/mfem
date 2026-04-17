#include "common.hpp"

extern "C" {

void CMFEM_PCG_SparseMatrixGSSmoother(const CMFEM_SparseMatrix *A,
                                      CMFEM_GSSmoother *M,
                                      const CMFEM_Vector *B,
                                      CMFEM_Vector *X,
                                      int print_iter,
                                      int max_iter,
                                      double rtol,
                                      double atol)
{
   mfem::PCG(*cmfem::As<const mfem::SparseMatrix>(A),
             *cmfem::As<mfem::GSSmoother>(M),
             cmfem::VectorRef(B),
             cmfem::VectorRef(X),
             print_iter,
             max_iter,
             static_cast<mfem::real_t>(rtol),
             static_cast<mfem::real_t>(atol));
}

void CMFEM_PCG_OperatorGSSmoother(const CMFEM_OperatorPtr *A,
                                  CMFEM_GSSmoother *M,
                                  const CMFEM_Vector *B,
                                  CMFEM_Vector *X,
                                  int print_iter,
                                  int max_iter,
                                  double rtol,
                                  double atol)
{
   mfem::PCG(*cmfem::OperatorPtrRef(A),
             *cmfem::As<mfem::GSSmoother>(M),
             cmfem::VectorRef(B),
             cmfem::VectorRef(X),
             print_iter,
             max_iter,
             static_cast<mfem::real_t>(rtol),
             static_cast<mfem::real_t>(atol));
}

void CMFEM_PCG_OperatorJacobiSmoother(const CMFEM_OperatorPtr *A,
                                      CMFEM_OperatorJacobiSmoother *M,
                                      const CMFEM_Vector *B,
                                      CMFEM_Vector *X,
                                      int print_iter,
                                      int max_iter,
                                      double rtol,
                                      double atol)
{
   mfem::PCG(*cmfem::OperatorPtrRef(A),
             *cmfem::As<mfem::OperatorJacobiSmoother>(M),
             cmfem::VectorRef(B),
             cmfem::VectorRef(X),
             print_iter,
             max_iter,
             static_cast<mfem::real_t>(rtol),
             static_cast<mfem::real_t>(atol));
}

void CMFEM_CG_Operator(const CMFEM_OperatorPtr *A,
                       const CMFEM_Vector *B,
                       CMFEM_Vector *X,
                       int print_iter,
                       int max_iter,
                       double rtol,
                       double atol)
{
   mfem::CG(*cmfem::OperatorPtrRef(A),
            cmfem::VectorRef(B),
            cmfem::VectorRef(X),
            print_iter,
            max_iter,
            static_cast<mfem::real_t>(rtol),
            static_cast<mfem::real_t>(atol));
}

int CMFEM_UsesTensorBasis(const CMFEM_FiniteElementSpace *fespace)
{
   return mfem::UsesTensorBasis(*cmfem::As<const mfem::FiniteElementSpace>(fespace));
}

} // extern "C"
