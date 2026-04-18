#include "common.hpp"

extern "C" {

   void CMFEM_PCGSmGs(const CMFEM_SparseMatrix *A,
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

   void CMFEM_PCGSmDs(const CMFEM_SparseMatrix *A,
                      CMFEM_DSmoother *M,
                      const CMFEM_Vector *B,
                      CMFEM_Vector *X,
                      int print_iter,
                      int max_iter,
                      double rtol,
                      double atol)
   {
      mfem::PCG(*cmfem::As<const mfem::SparseMatrix>(A),
                *cmfem::As<mfem::DSmoother>(M),
                cmfem::VectorRef(B),
                cmfem::VectorRef(X),
                print_iter,
                max_iter,
                static_cast<mfem::real_t>(rtol),
                static_cast<mfem::real_t>(atol));
   }

   void CMFEM_PCGBfOjs(const CMFEM_BilinearForm *A,
                       CMFEM_OperatorJacobiSmoother *M,
                       const CMFEM_Vector *B,
                       CMFEM_Vector *X,
                       int print_iter,
                       int max_iter,
                       double rtol,
                       double atol)
   {
      mfem::PCG(*cmfem::As<const mfem::BilinearForm>(A),
                *cmfem::As<mfem::OperatorJacobiSmoother>(M),
                cmfem::VectorRef(B),
                cmfem::VectorRef(X),
                print_iter,
                max_iter,
                static_cast<mfem::real_t>(rtol),
                static_cast<mfem::real_t>(atol));
   }

   void CMFEM_PCGBfDs(const CMFEM_BilinearForm *A,
                      CMFEM_DSmoother *M,
                      const CMFEM_Vector *B,
                      CMFEM_Vector *X,
                      int print_iter,
                      int max_iter,
                      double rtol,
                      double atol)
   {
      mfem::PCG(*cmfem::As<const mfem::BilinearForm>(A),
                *cmfem::As<mfem::DSmoother>(M),
                cmfem::VectorRef(B),
                cmfem::VectorRef(X),
                print_iter,
                max_iter,
                static_cast<mfem::real_t>(rtol),
                static_cast<mfem::real_t>(atol));
   }

   void CMFEM_PCGOpGs(const CMFEM_OperatorPtr *A,
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

   void CMFEM_PCGOpOjs(const CMFEM_OperatorPtr *A,
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

   void CMFEM_CGOp(const CMFEM_OperatorPtr *A,
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

   void CMFEM_GMRESOpGs(const CMFEM_OperatorPtr *A,
                        CMFEM_GSSmoother *M,
                        const CMFEM_Vector *B,
                        CMFEM_Vector *X,
                        int print_iter,
                        int max_iter,
                        int restart,
                        double rtol,
                        double atol)
   {
      mfem::GMRES(*cmfem::OperatorPtrRef(A),
                  *cmfem::As<mfem::GSSmoother>(M),
                  cmfem::VectorRef(B),
                  cmfem::VectorRef(X),
                  print_iter,
                  max_iter,
                  restart,
                  static_cast<mfem::real_t>(rtol),
                  static_cast<mfem::real_t>(atol));
   }

   void CMFEM_GMRESSmGs(const CMFEM_SparseMatrix *A,
                        CMFEM_GSSmoother *M,
                        const CMFEM_Vector *B,
                        CMFEM_Vector *X,
                        int print_iter,
                        int max_iter,
                        int restart,
                        double rtol,
                        double atol)
   {
      mfem::GMRES(*cmfem::As<const mfem::SparseMatrix>(A),
                  *cmfem::As<mfem::GSSmoother>(M),
                  cmfem::VectorRef(B),
                  cmfem::VectorRef(X),
                  print_iter,
                  max_iter,
                  restart,
                  static_cast<mfem::real_t>(rtol),
                  static_cast<mfem::real_t>(atol));
   }

   void CMFEM_PCGRAPOperatorBdp(const CMFEM_RAPOperator *A,
                                CMFEM_BlockDiagonalPreconditioner *M,
                                const CMFEM_Vector *B,
                                CMFEM_Vector *X,
                                int print_iter,
                                int max_iter,
                                double rtol,
                                double atol)
   {
      mfem::PCG(*cmfem::As<const mfem::RAPOperator>(A),
                *cmfem::As<mfem::BlockDiagonalPreconditioner>(M),
                cmfem::VectorRef(B),
                cmfem::VectorRef(X),
                print_iter,
                max_iter,
                static_cast<mfem::real_t>(rtol),
                static_cast<mfem::real_t>(atol));
   }

   void CMFEM_MINRESBopBdp(const CMFEM_BlockOperator *A,
                           CMFEM_BlockDiagonalPreconditioner *M,
                           const CMFEM_Vector *B,
                           CMFEM_Vector *X,
                           int print_iter,
                           int max_iter,
                           double rtol,
                           double atol)
   {
      const mfem::real_t rel_tol = static_cast<mfem::real_t>(rtol);
      const mfem::real_t abs_tol = static_cast<mfem::real_t>(atol);
      mfem::MINRES(*cmfem::As<const mfem::BlockOperator>(A),
                   *cmfem::As<mfem::BlockDiagonalPreconditioner>(M),
                   cmfem::VectorRef(B),
                   cmfem::VectorRef(X),
                   print_iter,
                   max_iter,
                   rel_tol * rel_tol,
                   abs_tol * abs_tol);
   }

   int CMFEM_UsesTensorBasis(const CMFEM_FiniteElementSpace *fespace)
   {
      return mfem::UsesTensorBasis(*cmfem::As<const mfem::FiniteElementSpace>
                                   (fespace));
   }

} // extern "C"
