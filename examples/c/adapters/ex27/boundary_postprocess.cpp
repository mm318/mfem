#include "common.hpp"

#include <cmath>

extern "C" {

   double CMFEM_IntegrateBoundaryConditionGfAi(const CMFEM_GridFunction
                                               *grid_function,
                                               const CMFEM_ArrayInt *marker,
                                               double alpha,
                                               double beta,
                                               double gamma,
                                               double *error)
   {
      double nrm = 0.0;
      double avg = 0.0;
      double local_error = 0.0;
      const bool alpha_is_zero = alpha == 0.0;
      const bool beta_is_zero = beta == 0.0;

      const mfem::GridFunction &x = *cmfem::As<const mfem::GridFunction>(
                                       grid_function);
      const mfem::FiniteElementSpace &fes = *x.FESpace();
      MFEM_VERIFY(fes.GetVDim() == 1,
                  "Boundary-condition integration expects a scalar field.");
      mfem::Mesh &mesh = *fes.GetMesh();
      const mfem::Array<int> &bdr = cmfem::ArrayIntRef(marker);

      mfem::Vector shape;
      mfem::Vector loc_dofs;
      mfem::Vector w_nor;
      mfem::DenseMatrix dshape;
      mfem::Array<int> dof_ids;

      for (int i = 0; i < mesh.GetNBE(); i++)
      {
         if (bdr[mesh.GetBdrAttribute(i) - 1] == 0)
         {
            continue;
         }

         mfem::FaceElementTransformations *ftr = mesh.GetBdrFaceTransformations(i);
         if (ftr == nullptr)
         {
            continue;
         }

         const mfem::FiniteElement &fe = *fes.GetFE(ftr->Elem1No);
         MFEM_VERIFY(fe.GetMapType() == mfem::FiniteElement::VALUE,
                     "Boundary-condition integration expects VALUE finite elements.");
         const int int_order = 2 * fe.GetOrder() + 3;
         const mfem::IntegrationRule &ir = mfem::IntRules.Get(ftr->FaceGeom,
                                                              int_order);

         fes.GetElementDofs(ftr->Elem1No, dof_ids);
         x.GetSubVector(dof_ids, loc_dofs);
         if (!alpha_is_zero)
         {
            const int sdim = ftr->Face->GetSpaceDim();
            w_nor.SetSize(sdim);
            dshape.SetSize(fe.GetDof(), sdim);
         }
         if (!beta_is_zero)
         {
            shape.SetSize(fe.GetDof());
         }

         for (int j = 0; j < ir.GetNPoints(); j++)
         {
            const mfem::IntegrationPoint &ip = ir.IntPoint(j);
            mfem::IntegrationPoint eip;
            ftr->Loc1.Transform(ip, eip);
            ftr->Face->SetIntPoint(&ip);
            const double face_weight = ftr->Face->Weight();
            double val = 0.0;

            if (!alpha_is_zero)
            {
               ftr->Elem1->SetIntPoint(&eip);
               fe.CalcPhysDShape(*ftr->Elem1, dshape);
               mfem::CalcOrtho(ftr->Face->Jacobian(), w_nor);
               val += alpha * dshape.InnerProduct(w_nor, loc_dofs) / face_weight;
            }
            if (!beta_is_zero)
            {
               fe.CalcShape(eip, shape);
               val += beta * (shape * loc_dofs);
            }

            nrm += ip.weight * face_weight;
            avg += val * ip.weight * face_weight;

            val -= gamma;
            local_error += (val * val) * ip.weight * face_weight;
         }
      }

      if (std::abs(nrm) > 0.0)
      {
         local_error /= nrm;
         avg /= nrm;
      }

      local_error = local_error >= 0.0 ?
                    std::sqrt(local_error) : -std::sqrt(-local_error);
      if (error)
      {
         *error = local_error;
      }
      return avg;
   }

} // extern "C"
