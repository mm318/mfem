#include "common.hpp"

extern "C" {

   double CMFEM_ComputeLpNormFcMeshOrder(
      double p,
      const CMFEM_FunctionCoefficient *coefficient,
      const CMFEM_Mesh *mesh,
      int quadrature_order)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
      for (int i = 0; i < mfem::Geometry::NumGeom; i++)
      {
         irs[i] = &mfem::IntRules.Get(i, quadrature_order);
      }
      return mfem::ComputeLpNorm(static_cast<mfem::real_t>(p),
                                 coef,
                                 const_cast<mfem::Mesh &>(
                                    *cmfem::As<const mfem::Mesh>(mesh)),
                                 irs);
   }

   double CMFEM_ComputeLpNormVfcMeshOrder(
      double p,
      const CMFEM_VectorFunctionCoefficient *coefficient,
      const CMFEM_Mesh *mesh,
      int quadrature_order)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
      for (int i = 0; i < mfem::Geometry::NumGeom; i++)
      {
         irs[i] = &mfem::IntRules.Get(i, quadrature_order);
      }
      return mfem::ComputeLpNorm(static_cast<mfem::real_t>(p),
                                 coef,
                                 const_cast<mfem::Mesh &>(
                                    *cmfem::As<const mfem::Mesh>(mesh)),
                                 irs);
   }

} // extern "C"
