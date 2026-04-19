#include "common.hpp"

#include <complex>
#include <fstream>

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ComplexGridFunction, mfem::ComplexGridFunction);

const mfem::Array<int> *ToElementMarker(const CMFEM_ArrayInt *element_marker)
{
   return element_marker ? &cmfem::ArrayIntRef(element_marker) : nullptr;
}

const mfem::IntegrationRule * const *BuildIntegrationRules(
   int quadrature_order,
   const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom])
{
   for (int i = 0; i < mfem::Geometry::NumGeom; i++)
   {
      irs[i] = &mfem::IntRules.Get(i, quadrature_order);
   }
   return irs;
}

} // namespace

extern "C" {

   CMFEM_ComplexGridFunction *CMFEM_ComplexGridFunction_New(
      CMFEM_FiniteElementSpace *fespace)
   {
      return reinterpret_cast<CMFEM_ComplexGridFunction *>(
                new mfem::ComplexGridFunction(
                   cmfem::As<mfem::FiniteElementSpace>(fespace)));
   }

   void CMFEM_ComplexGridFunction_Delete(
      CMFEM_ComplexGridFunction *grid_function)
   {
      delete cmfem::As<mfem::ComplexGridFunction>(grid_function);
   }

   void CMFEM_ComplexGridFunction_Assign(
      CMFEM_ComplexGridFunction *grid_function,
      double real_value,
      double imag_value)
   {
      *cmfem::As<mfem::ComplexGridFunction>(grid_function) =
         std::complex<mfem::real_t>(static_cast<mfem::real_t>(real_value),
                                    static_cast<mfem::real_t>(imag_value));
   }

   void CMFEM_ComplexGridFunction_ProjectCoefficientFcFc(
      CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_FunctionCoefficient *real_coefficient,
      const CMFEM_FunctionCoefficient *imag_coefficient)
   {
      auto &real_coef = const_cast<mfem::FunctionCoefficient &>(
                           *cmfem::As<const mfem::FunctionCoefficient>(real_coefficient));
      auto &imag_coef = const_cast<mfem::FunctionCoefficient &>(
                           *cmfem::As<const mfem::FunctionCoefficient>(imag_coefficient));
      cmfem::As<mfem::ComplexGridFunction>(grid_function)->ProjectCoefficient(
         real_coef,
         imag_coef);
   }

   void CMFEM_ComplexGridFunction_ProjectCoefficientVfcVfc(
      CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *real_coefficient,
      const CMFEM_VectorFunctionCoefficient *imag_coefficient)
   {
      auto &real_coef = const_cast<mfem::VectorFunctionCoefficient &>(
                           *cmfem::As<const mfem::VectorFunctionCoefficient>(real_coefficient));
      auto &imag_coef = const_cast<mfem::VectorFunctionCoefficient &>(
                           *cmfem::As<const mfem::VectorFunctionCoefficient>(imag_coefficient));
      cmfem::As<mfem::ComplexGridFunction>(grid_function)->ProjectCoefficient(
         real_coef,
         imag_coef);
   }

   void CMFEM_ComplexGridFunction_ProjectBdrCoefficientFcFcAi(
      CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_FunctionCoefficient *real_coefficient,
      const CMFEM_FunctionCoefficient *imag_coefficient,
      const CMFEM_ArrayInt *attributes)
   {
      auto &real_coef = const_cast<mfem::FunctionCoefficient &>(
                           *cmfem::As<const mfem::FunctionCoefficient>(real_coefficient));
      auto &imag_coef = const_cast<mfem::FunctionCoefficient &>(
                           *cmfem::As<const mfem::FunctionCoefficient>(imag_coefficient));
      cmfem::As<mfem::ComplexGridFunction>(grid_function)->ProjectBdrCoefficient(
         real_coef,
         imag_coef,
         const_cast<mfem::Array<int> &>(cmfem::ArrayIntRef(attributes)));
   }

   void CMFEM_ComplexGridFunction_ProjectBdrCoefficientTangentVfcVfcAi(
      CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *real_coefficient,
      const CMFEM_VectorFunctionCoefficient *imag_coefficient,
      const CMFEM_ArrayInt *attributes)
   {
      auto &real_coef = const_cast<mfem::VectorFunctionCoefficient &>(
                           *cmfem::As<const mfem::VectorFunctionCoefficient>(real_coefficient));
      auto &imag_coef = const_cast<mfem::VectorFunctionCoefficient &>(
                           *cmfem::As<const mfem::VectorFunctionCoefficient>(imag_coefficient));
      cmfem::As<mfem::ComplexGridFunction>(grid_function)
      ->ProjectBdrCoefficientTangent(real_coef, imag_coef,
                                     const_cast<mfem::Array<int> &>(
                                        cmfem::ArrayIntRef(attributes)));
   }

   void CMFEM_ComplexGridFunction_ProjectBdrCoefficientNormalVfcVfcAi(
      CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *real_coefficient,
      const CMFEM_VectorFunctionCoefficient *imag_coefficient,
      const CMFEM_ArrayInt *attributes)
   {
      auto &real_coef = const_cast<mfem::VectorFunctionCoefficient &>(
                           *cmfem::As<const mfem::VectorFunctionCoefficient>(real_coefficient));
      auto &imag_coef = const_cast<mfem::VectorFunctionCoefficient &>(
                           *cmfem::As<const mfem::VectorFunctionCoefficient>(imag_coefficient));
      cmfem::As<mfem::ComplexGridFunction>(grid_function)
      ->ProjectBdrCoefficientNormal(real_coef, imag_coef,
                                    const_cast<mfem::Array<int> &>(
                                       cmfem::ArrayIntRef(attributes)));
   }

   double CMFEM_ComplexGridFunction_ComputeL2ErrorRealFc(
      const CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_FunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::ComplexGridFunction>(grid_function)
             ->real().ComputeL2Error(coef);
   }

   double CMFEM_ComplexGridFunction_ComputeL2ErrorImagFc(
      const CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_FunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::ComplexGridFunction>(grid_function)
             ->imag().ComputeL2Error(coef);
   }

   double CMFEM_ComplexGridFunction_ComputeL2ErrorRealVfc(
      const CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::ComplexGridFunction>(grid_function)
             ->real().ComputeL2Error(coef);
   }

   double CMFEM_ComplexGridFunction_ComputeL2ErrorImagVfc(
      const CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      return cmfem::As<const mfem::ComplexGridFunction>(grid_function)
             ->imag().ComputeL2Error(coef);
   }

   double CMFEM_ComplexGridFunction_ComputeL2ErrorRealVfcAiQo(
      const CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient,
      const CMFEM_ArrayInt *element_marker,
      int quadrature_order)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
      auto *nonconst_irs =
         const_cast<const mfem::IntegrationRule **>(
            BuildIntegrationRules(quadrature_order, irs));
      return cmfem::As<const mfem::ComplexGridFunction>(grid_function)
             ->real().ComputeL2Error(coef,
                                     nonconst_irs,
                                     ToElementMarker(element_marker));
   }

   double CMFEM_ComplexGridFunction_ComputeL2ErrorImagVfcAiQo(
      const CMFEM_ComplexGridFunction *grid_function,
      const CMFEM_VectorFunctionCoefficient *coefficient,
      const CMFEM_ArrayInt *element_marker,
      int quadrature_order)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
      auto *nonconst_irs =
         const_cast<const mfem::IntegrationRule **>(
            BuildIntegrationRules(quadrature_order, irs));
      return cmfem::As<const mfem::ComplexGridFunction>(grid_function)
             ->imag().ComputeL2Error(coef,
                                     nonconst_irs,
                                     ToElementMarker(element_marker));
   }

   void CMFEM_ComplexGridFunction_CopyRealToGf(
      const CMFEM_ComplexGridFunction *grid_function,
      CMFEM_GridFunction *real_grid_function)
   {
      *cmfem::As<mfem::GridFunction>(real_grid_function) =
         cmfem::As<const mfem::ComplexGridFunction>(grid_function)->real();
   }

   void CMFEM_ComplexGridFunction_CopyImagToGf(
      const CMFEM_ComplexGridFunction *grid_function,
      CMFEM_GridFunction *imag_grid_function)
   {
      *cmfem::As<mfem::GridFunction>(imag_grid_function) =
         cmfem::As<const mfem::ComplexGridFunction>(grid_function)->imag();
   }

   void CMFEM_ComplexGridFunction_Save(
      const CMFEM_ComplexGridFunction *grid_function,
      const char *path,
      int precision)
   {
      std::ofstream out(path);
      out.precision(precision);
      cmfem::As<const mfem::ComplexGridFunction>(grid_function)->Save(out);
   }

   void CMFEM_ComplexGridFunction_SaveReal(
      const CMFEM_ComplexGridFunction *grid_function,
      const char *path,
      int precision)
   {
      std::ofstream out(path);
      out.precision(precision);
      cmfem::As<const mfem::ComplexGridFunction>(grid_function)->real().Save(out);
   }

   void CMFEM_ComplexGridFunction_SaveImag(
      const CMFEM_ComplexGridFunction *grid_function,
      const char *path,
      int precision)
   {
      std::ofstream out(path);
      out.precision(precision);
      cmfem::As<const mfem::ComplexGridFunction>(grid_function)->imag().Save(out);
   }

} // extern "C"
