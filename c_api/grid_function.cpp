#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_GridFunction, mfem::GridFunction);

} // namespace

extern "C" {

   CMFEM_GridFunction *CMFEM_GridFunction_New(CMFEM_FiniteElementSpace *fespace)
   {
      return reinterpret_cast<CMFEM_GridFunction *>(
                new mfem::GridFunction(cmfem::As<mfem::FiniteElementSpace>(fespace)));
   }

   void CMFEM_GridFunction_Delete(CMFEM_GridFunction *grid_function)
   {
      delete cmfem::As<mfem::GridFunction>(grid_function);
   }

   void CMFEM_GridFunction_Assign(CMFEM_GridFunction *grid_function, double value)
   {
      *cmfem::As<mfem::GridFunction>(grid_function) = static_cast<mfem::real_t>
                                                      (value);
   }

   void CMFEM_GridFunction_ProjectVectorFunctionCoefficient(
      CMFEM_GridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      cmfem::As<mfem::GridFunction>(grid_function)->ProjectCoefficient(coef);
   }

   double CMFEM_GridFunction_ComputeL2ErrorFunctionCoefficient(
      const CMFEM_GridFunction *grid_function,
      const CMFEM_FunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeL2Error(coef);
   }

   double CMFEM_GridFunction_ComputeL2ErrorVectorFunctionCoefficient(
      const CMFEM_GridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeL2Error(coef);
   }

   void CMFEM_GridFunction_Add(CMFEM_GridFunction *grid_function,
                               const CMFEM_GridFunction *other)
   {
      *cmfem::As<mfem::GridFunction>(grid_function) +=
         *cmfem::As<const mfem::GridFunction>(other);
   }

   void CMFEM_GridFunction_Scale(CMFEM_GridFunction *grid_function, double scale)
   {
      *cmfem::As<mfem::GridFunction>(grid_function) *= static_cast<mfem::real_t>
                                                       (scale);
   }

   void CMFEM_GridFunction_Save(const CMFEM_GridFunction *grid_function,
                                const char *path,
                                int precision)
   {
      cmfem::As<const mfem::GridFunction>(grid_function)->Save(path, precision);
   }

} // extern "C"
