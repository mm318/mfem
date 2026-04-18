#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_SparseMatrix, mfem::SparseMatrix);

} // namespace

extern "C" {

   CMFEM_SparseMatrix CMFEM_SparseMatrix_Construct(void)
   {
      alignas(mfem::SparseMatrix) CMFEM_SparseMatrix matrix;
      new (cmfem::As<mfem::SparseMatrix>(&matrix)) mfem::SparseMatrix();
      return matrix;
   }

   CMFEM_SparseMatrix *CMFEM_SparseMatrix_New(void)
   {
      return reinterpret_cast<CMFEM_SparseMatrix *>(new mfem::SparseMatrix());
   }

   CMFEM_SparseMatrix CMFEM_SparseMatrix_Copy(const CMFEM_SparseMatrix *matrix)
   {
      alignas(mfem::SparseMatrix) CMFEM_SparseMatrix copy;
      new (cmfem::As<mfem::SparseMatrix>(&copy)) mfem::SparseMatrix(
         *cmfem::As<const mfem::SparseMatrix>(matrix));
      return copy;
   }

   CMFEM_SparseMatrix *CMFEM_SparseMatrix_NewCopy(const CMFEM_SparseMatrix *matrix)
   {
      return reinterpret_cast<CMFEM_SparseMatrix *>(
                new mfem::SparseMatrix(*cmfem::As<const mfem::SparseMatrix>(matrix)));
   }

   void CMFEM_SparseMatrix_Destroy(CMFEM_SparseMatrix *matrix)
   {
      cmfem::As<mfem::SparseMatrix>(matrix)->~SparseMatrix();
   }

   void CMFEM_SparseMatrix_Delete(CMFEM_SparseMatrix *matrix)
   {
      delete cmfem::As<mfem::SparseMatrix>(matrix);
   }

   int CMFEM_SparseMatrix_Height(const CMFEM_SparseMatrix *matrix)
   {
      return cmfem::As<const mfem::SparseMatrix>(matrix)->Height();
   }

   void CMFEM_SparseMatrix_Mult(const CMFEM_SparseMatrix *matrix,
                                const CMFEM_Vector *x,
                                CMFEM_Vector *y)
   {
      cmfem::As<const mfem::SparseMatrix>(matrix)->Mult(cmfem::VectorRef(x),
                                                        cmfem::VectorRef(y));
   }

   double CMFEM_SparseMatrix_InnerProduct(const CMFEM_SparseMatrix *matrix,
                                          const CMFEM_Vector *x,
                                          const CMFEM_Vector *y)
   {
      return cmfem::As<const mfem::SparseMatrix>(matrix)->InnerProduct(
                cmfem::VectorRef(x),
                cmfem::VectorRef(y));
   }

   CMFEM_SparseMatrix *CMFEM_RAPSmSmSm(const CMFEM_SparseMatrix *Rt,
                                       const CMFEM_SparseMatrix *A,
                                       const CMFEM_SparseMatrix *P)
   {
      return reinterpret_cast<CMFEM_SparseMatrix *>(
                mfem::RAP(*cmfem::As<const mfem::SparseMatrix>(Rt),
                          *cmfem::As<const mfem::SparseMatrix>(A),
                          *cmfem::As<const mfem::SparseMatrix>(P)));
   }

   CMFEM_SparseMatrix *CMFEM_AddSmSm(double a,
                                     const CMFEM_SparseMatrix *A,
                                     double b,
                                     const CMFEM_SparseMatrix *B)
   {
      return reinterpret_cast<CMFEM_SparseMatrix *>(mfem::Add(
                                                       static_cast<mfem::real_t>(
                                                          a),
                                                       *cmfem::As<const mfem::SparseMatrix>(
                                                          A),
                                                       static_cast<mfem::real_t>(
                                                          b),
                                                       *cmfem::As<const mfem::SparseMatrix>(
                                                          B)));
   }

} // extern "C"
