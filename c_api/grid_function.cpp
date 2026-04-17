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

   void CMFEM_GridFunction_Set(CMFEM_GridFunction *grid_function,
                               int index,
                               double value)
   {
      (*cmfem::As<mfem::GridFunction>(grid_function))(index) =
         static_cast<mfem::real_t>(value);
   }

   void CMFEM_GridFunction_ProjectCoefficientFc(
      CMFEM_GridFunction *grid_function,
      const CMFEM_FunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      cmfem::As<mfem::GridFunction>(grid_function)->ProjectCoefficient(coef);
   }

   void CMFEM_GridFunction_ProjectBdrCoefficientCc(
      CMFEM_GridFunction *grid_function,
      const CMFEM_ConstantCoefficient *coefficient,
      const CMFEM_ArrayInt *attributes)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::GridFunction>(grid_function)->ProjectBdrCoefficient(
         coef,
         cmfem::ArrayIntRef(attributes));
   }

   void CMFEM_GridFunction_ProjectBdrCoefficientFc(
      CMFEM_GridFunction *grid_function,
      const CMFEM_FunctionCoefficient *coefficient,
      const CMFEM_ArrayInt *attributes)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      cmfem::As<mfem::GridFunction>(grid_function)->ProjectBdrCoefficient(
         coef,
         cmfem::ArrayIntRef(attributes));
   }

   void CMFEM_GridFunction_ProjectBdrCoefficientVcc(
      CMFEM_GridFunction *grid_function,
      const CMFEM_VectorConstantCoefficient *coefficient,
      const CMFEM_ArrayInt *attributes)
   {
      auto &coef = const_cast<mfem::VectorConstantCoefficient &>(
                      *cmfem::As<const mfem::VectorConstantCoefficient>(coefficient));
      cmfem::As<mfem::GridFunction>(grid_function)->ProjectBdrCoefficient(
         coef,
         cmfem::ArrayIntRef(attributes));
   }

   void CMFEM_GridFunction_ProjectCoefficientVfc(
      CMFEM_GridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      cmfem::As<mfem::GridFunction>(grid_function)->ProjectCoefficient(coef);
   }

   double CMFEM_GridFunction_ComputeL2ErrorFc(
      const CMFEM_GridFunction *grid_function,
      const CMFEM_FunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeL2Error(coef);
   }

   double CMFEM_GridFunction_ComputeL2ErrorVfc(
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

   void CMFEM_GridFunction_GetTrueDofs(const CMFEM_GridFunction *grid_function,
                                       CMFEM_Vector *true_dofs)
   {
      cmfem::As<const mfem::GridFunction>(grid_function)->GetTrueDofs(
         cmfem::VectorRef(true_dofs));
   }

   void CMFEM_GridFunction_SetFromTrueDofs(CMFEM_GridFunction *grid_function,
                                           const CMFEM_Vector *true_dofs)
   {
      cmfem::As<mfem::GridFunction>(grid_function)->SetFromTrueDofs(
         cmfem::VectorRef(true_dofs));
   }

   void CMFEM_GridFunction_Update(CMFEM_GridFunction *grid_function)
   {
      cmfem::As<mfem::GridFunction>(grid_function)->Update();
   }

   void CMFEM_GridFunction_Save(const CMFEM_GridFunction *grid_function,
                                const char *path,
                                int precision)
   {
      cmfem::As<const mfem::GridFunction>(grid_function)->Save(path, precision);
   }

} // extern "C"
