#include "common.hpp"

#include <vector>

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_MatrixFunctionCoefficient,
                  mfem::MatrixFunctionCoefficient);

} // namespace

extern "C" {

   CMFEM_MatrixFunctionCoefficient *CMFEM_MatrixFunctionCoefficient_New(
      int dim,
      CMFEM_MatrixFunctionCallback callback,
      void *context)
   {
      auto function = [callback, context, dim](const mfem::Vector &x,
                                               mfem::DenseMatrix &value)
      {
         std::vector<double> entries(static_cast<size_t>(dim * dim));
         callback(reinterpret_cast<const CMFEM_Vector *>(&x),
                  entries.data(),
                  dim,
                  context);
         value.SetSize(dim);
         for (int row = 0; row < dim; row++)
         {
            for (int col = 0; col < dim; col++)
            {
               value(row, col) = static_cast<mfem::real_t>(
                                    entries[static_cast<size_t>(row * dim + col)]);
            }
         }
      };
      return reinterpret_cast<CMFEM_MatrixFunctionCoefficient *>(
                new mfem::MatrixFunctionCoefficient(dim, std::move(function)));
   }

   void CMFEM_MatrixFunctionCoefficient_Delete(
      CMFEM_MatrixFunctionCoefficient *coefficient)
   {
      delete cmfem::As<mfem::MatrixFunctionCoefficient>(coefficient);
   }

} // extern "C"
