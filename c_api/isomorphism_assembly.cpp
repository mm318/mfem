#include "common.hpp"

namespace
{

class IsomorphismCoefficient : public mfem::VectorCoefficient
{
private:
   mfem::GridFunction &psi;

public:
   IsomorphismCoefficient(int vdim, mfem::GridFunction &psi_)
      : mfem::VectorCoefficient(vdim), psi(psi_) { }

   using mfem::VectorCoefficient::Eval;

   void Eval(mfem::Vector &value,
             mfem::ElementTransformation &transform,
             const mfem::IntegrationPoint &ip) override
   {
      mfem::Vector psi_vals(vdim);
      psi.GetVectorValue(transform, ip, psi_vals);
      const mfem::real_t norm = psi_vals.Norml2();
      const mfem::real_t phi = 1.0 / std::sqrt(1.0 + norm * norm);

      value = psi_vals;
      value *= phi;
   }
};

class DIsomorphismCoefficient : public mfem::MatrixCoefficient
{
private:
   mfem::GridFunction &psi;
   mfem::real_t eps;

public:
   DIsomorphismCoefficient(int height, mfem::GridFunction &psi_, mfem::real_t eps_)
      : mfem::MatrixCoefficient(height), psi(psi_), eps(eps_) { }

   void Eval(mfem::DenseMatrix &value,
             mfem::ElementTransformation &transform,
             const mfem::IntegrationPoint &ip) override
   {
      mfem::Vector psi_vals(height);
      psi.GetVectorValue(transform, ip, psi_vals);
      const mfem::real_t norm = psi_vals.Norml2();
      const mfem::real_t phi = 1.0 / std::sqrt(1.0 + norm * norm);

      value.SetSize(height);
      value = 0.0;
      for (int i = 0; i < height; i++)
      {
         value(i, i) = phi + eps;
         for (int j = 0; j < height; j++)
         {
            value(i, j) -= psi_vals(i) * psi_vals(j) * std::pow(phi, 3);
         }
      }
   }
};

} // namespace

extern "C" {

   void CMFEM_AssembleIsomorphismLinearFormRtGf(
      const CMFEM_FiniteElementSpace *fespace,
      const CMFEM_GridFunction *psi,
      CMFEM_Vector *out)
   {
      auto &space = *cmfem::As<const mfem::FiniteElementSpace>(fespace);
      auto &psi_ref = const_cast<mfem::GridFunction &>(
                         *cmfem::As<const mfem::GridFunction>(psi));
      IsomorphismCoefficient coeff(space.GetMesh()->SpaceDimension(), psi_ref);
      mfem::LinearForm linear_form(const_cast<mfem::FiniteElementSpace *>(&space));
      linear_form.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(coeff));
      linear_form.Assemble();
      cmfem::VectorRef(out) = linear_form;
   }

   CMFEM_SparseMatrix *CMFEM_AssembleDIsomorphismMassMatrixRtGf(
      const CMFEM_FiniteElementSpace *fespace,
      const CMFEM_GridFunction *psi,
      double eps)
   {
      auto &space = *cmfem::As<const mfem::FiniteElementSpace>(fespace);
      auto &psi_ref = const_cast<mfem::GridFunction &>(
                         *cmfem::As<const mfem::GridFunction>(psi));
      DIsomorphismCoefficient coeff(space.GetMesh()->SpaceDimension(),
                                    psi_ref,
                                    static_cast<mfem::real_t>(eps));
      mfem::BilinearForm bilinear_form(const_cast<mfem::FiniteElementSpace *>
                                       (&space));
      bilinear_form.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(coeff));
      bilinear_form.Assemble(false);
      bilinear_form.Finalize(false);
      return reinterpret_cast<CMFEM_SparseMatrix *>(
                new mfem::SparseMatrix(bilinear_form.SpMat()));
   }

} // extern "C"
