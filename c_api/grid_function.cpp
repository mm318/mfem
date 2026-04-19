#include "common.hpp"

#include <fstream>

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

   double CMFEM_GridFunction_Get(const CMFEM_GridFunction *grid_function,
                                 int index)
   {
      return (*cmfem::As<const mfem::GridFunction>(grid_function))(index);
   }

   void CMFEM_GridFunction_Set(CMFEM_GridFunction *grid_function,
                               int index,
                               double value)
   {
      (*cmfem::As<mfem::GridFunction>(grid_function))(index) =
         static_cast<mfem::real_t>(value);
   }

   int CMFEM_GridFunction_Size(const CMFEM_GridFunction *grid_function)
   {
      return cmfem::As<const mfem::GridFunction>(grid_function)->Size();
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

   void CMFEM_GridFunction_ComputeFluxDi(
      const CMFEM_GridFunction *grid_function,
      CMFEM_DiffusionIntegrator *integrator,
      CMFEM_GridFunction *flux)
   {
      const_cast<mfem::GridFunction *>(
         cmfem::As<const mfem::GridFunction>(grid_function))->ComputeFlux(
            *cmfem::As<mfem::DiffusionIntegrator>(integrator),
            *cmfem::As<mfem::GridFunction>(flux));
   }

   void CMFEM_GridFunction_SetTrueVector(CMFEM_GridFunction *grid_function)
   {
      cmfem::As<mfem::GridFunction>(grid_function)->SetTrueVector();
   }

   void CMFEM_GridFunction_SetFromTrueVector(CMFEM_GridFunction *grid_function)
   {
      cmfem::As<mfem::GridFunction>(grid_function)->SetFromTrueVector();
   }

   double CMFEM_GridFunction_ComputeL2ErrorFc(
      const CMFEM_GridFunction *grid_function,
      const CMFEM_FunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeL2Error(coef);
   }

   double CMFEM_GridFunction_ComputeL2ErrorFcOrder(
      const CMFEM_GridFunction *grid_function,
      const CMFEM_FunctionCoefficient *coefficient,
      int quadrature_order)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
      for (int i = 0; i < mfem::Geometry::NumGeom; i++)
      {
         irs[i] = &mfem::IntRules.Get(i, quadrature_order);
      }
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeL2Error(
                coef,
                irs);
   }

   double CMFEM_GridFunction_ComputeL2ErrorVfc(
      const CMFEM_GridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeL2Error(coef);
   }

   double CMFEM_GridFunction_ComputeL2ErrorVfcOrder(
      const CMFEM_GridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient,
      int quadrature_order)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
      for (int i = 0; i < mfem::Geometry::NumGeom; i++)
      {
         irs[i] = &mfem::IntRules.Get(i, quadrature_order);
      }
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeL2Error(
                coef,
                irs);
   }

   double CMFEM_GridFunction_ComputeHCurlErrorVfcVfc(
      const CMFEM_GridFunction *grid_function,
      CMFEM_VectorFunctionCoefficient *exact_solution,
      CMFEM_VectorFunctionCoefficient *exact_curl)
   {
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeHCurlError(
                cmfem::As<mfem::VectorFunctionCoefficient>(exact_solution),
                cmfem::As<mfem::VectorFunctionCoefficient>(exact_curl));
   }

   double CMFEM_GridFunction_ComputeH1ErrorFcVfc(
      const CMFEM_GridFunction *grid_function,
      CMFEM_FunctionCoefficient *exact_solution,
      CMFEM_VectorFunctionCoefficient *exact_gradient)
   {
      return cmfem::As<const mfem::GridFunction>(grid_function)->ComputeH1Error(
                cmfem::As<mfem::FunctionCoefficient>(exact_solution),
                cmfem::As<mfem::VectorFunctionCoefficient>(exact_gradient));
   }

   void CMFEM_GridFunction_CopyToVec(const CMFEM_GridFunction *grid_function,
                                     CMFEM_Vector *vector)
   {
      cmfem::VectorRef(vector) = *cmfem::As<const mfem::GridFunction>(grid_function);
   }

   void CMFEM_GridFunction_SetFromVec(CMFEM_GridFunction *grid_function,
                                      const CMFEM_Vector *vector)
   {
      *cmfem::As<mfem::GridFunction>(grid_function) = cmfem::VectorRef(vector);
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
      std::ofstream out(path);
      out.precision(precision);
      cmfem::As<const mfem::GridFunction>(grid_function)->Save(out);
   }

} // extern "C"
