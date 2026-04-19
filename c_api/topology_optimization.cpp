#include "common.hpp"

#include "examples/cpp/ex37.hpp"

#include <cmath>
#include <memory>

namespace
{

mfem::real_t Ex37ProjectVolume(mfem::GridFunction &psi,
                               mfem::real_t target_volume,
                               mfem::real_t tol,
                               int max_its)
{
   mfem::MappedGridFunctionCoefficient sigmoid_psi(&psi, mfem::sigmoid);
   mfem::MappedGridFunctionCoefficient der_sigmoid_psi(&psi, mfem::der_sigmoid);

   mfem::LinearForm int_sigmoid_psi(psi.FESpace());
   int_sigmoid_psi.AddDomainIntegrator(new mfem::DomainLFIntegrator(sigmoid_psi));
   mfem::LinearForm int_der_sigmoid_psi(psi.FESpace());
   int_der_sigmoid_psi.AddDomainIntegrator(
      new mfem::DomainLFIntegrator(der_sigmoid_psi));

   bool done = false;
   for (int k = 0; k < max_its; k++)
   {
      int_sigmoid_psi.Assemble();
      const mfem::real_t f = int_sigmoid_psi.Sum() - target_volume;

      int_der_sigmoid_psi.Assemble();
      const mfem::real_t df = int_der_sigmoid_psi.Sum();

      const mfem::real_t dc = -f / df;
      psi += dc;
      if (std::abs(dc) < tol)
      {
         done = true;
         break;
      }
   }
   if (!done)
   {
      mfem::mfem_warning("Projection reached maximum iteration without converging. "
                         "Result may not be accurate.");
   }
   int_sigmoid_psi.Assemble();
   return int_sigmoid_psi.Sum();
}

class Ex37FilterSolverContext
{
private:
   mfem::ConstantCoefficient epsilon_sq_;
   mfem::ConstantCoefficient one_;
   mfem::Array<int> ess_bdr_;
   mfem::DiffusionSolver solver_;

public:
   Ex37FilterSolverContext(mfem::Mesh *mesh,
                           int order,
                           mfem::real_t epsilon_sq,
                           const mfem::Array<int> *ess_bdr)
      : epsilon_sq_(epsilon_sq),
        one_(1.0)
   {
      solver_.SetMesh(mesh);
      solver_.SetOrder(order);
      solver_.SetDiffusionCoefficient(&epsilon_sq_);
      solver_.SetMassCoefficient(&one_);
      if (ess_bdr != nullptr)
      {
         ess_bdr_ = *ess_bdr;
      }
      solver_.SetEssentialBoundary(ess_bdr_);
      solver_.SetupFEM();
   }

   void SolveSigmoid(const mfem::GridFunction &psi, mfem::GridFunction &rho_filter)
   {
      mfem::MappedGridFunctionCoefficient rho(&psi, mfem::sigmoid);
      solver_.SetRHSCoefficient(&rho);
      solver_.Solve();
      rho_filter = *solver_.GetFEMSolution();
   }

   void SolveStrainEnergy(mfem::real_t lambda,
                          mfem::real_t mu,
                          const mfem::GridFunction &u,
                          const mfem::GridFunction &rho_filter,
                          mfem::real_t rho_min,
                          mfem::GridFunction &w_filter)
   {
      mfem::ConstantCoefficient lambda_cf(lambda);
      mfem::ConstantCoefficient mu_cf(mu);
      mfem::StrainEnergyDensityCoefficient rhs_cf(&lambda_cf,
                                                  &mu_cf,
                                                  const_cast<mfem::GridFunction *>(&u),
                                                  const_cast<mfem::GridFunction *>(&rho_filter),
                                                  rho_min);
      solver_.SetRHSCoefficient(&rhs_cf);
      solver_.Solve();
      w_filter = *solver_.GetFEMSolution();
   }
};

class Ex37ElasticitySolverContext
{
private:
   mfem::Vector center_;
   mfem::Vector force_;
   mfem::VolumeForceCoefficient rhs_cf_;
   mfem::Array<int> ess_bdr_;
   mfem::LinearElasticitySolver solver_;

public:
   Ex37ElasticitySolverContext(mfem::Mesh *mesh,
                               int order,
                               const mfem::Array<int> *ess_bdr,
                               mfem::real_t radius,
                               const double *center_xy,
                               const double *force_xy)
      : center_(2),
        force_(2),
        rhs_cf_(radius, center_, force_)
   {
      center_(0) = static_cast<mfem::real_t>(center_xy[0]);
      center_(1) = static_cast<mfem::real_t>(center_xy[1]);
      force_(0) = static_cast<mfem::real_t>(force_xy[0]);
      force_(1) = static_cast<mfem::real_t>(force_xy[1]);

      solver_.SetMesh(mesh);
      solver_.SetOrder(order);
      solver_.SetRHSCoefficient(&rhs_cf_);
      if (ess_bdr != nullptr)
      {
         ess_bdr_ = *ess_bdr;
      }
      solver_.SetEssentialBoundary(ess_bdr_);
      solver_.SetupFEM();
   }

   void Solve(mfem::real_t lambda,
              mfem::real_t mu,
              const mfem::GridFunction &rho_filter,
              mfem::real_t rho_min,
              mfem::GridFunction &u)
   {
      mfem::ConstantCoefficient lambda_cf(lambda);
      mfem::ConstantCoefficient mu_cf(mu);
      mfem::SIMPInterpolationCoefficient simp_cf(
         const_cast<mfem::GridFunction *>(&rho_filter), rho_min, 1.0);
      mfem::ProductCoefficient lambda_simp_cf(lambda_cf, simp_cf);
      mfem::ProductCoefficient mu_simp_cf(mu_cf, simp_cf);
      solver_.SetLameCoefficients(&lambda_simp_cf, &mu_simp_cf);
      solver_.Solve();
      u = *solver_.GetFEMSolution();
   }

   mfem::real_t Compliance(const mfem::GridFunction &u)
   {
      return (*solver_.GetLinearForm())(u);
   }
};

} // namespace

extern "C" {

   CMFEM_Mesh *CMFEM_Ex37_NewBeamMesh(void)
   {
      auto *mesh = new mfem::Mesh(
         mfem::Mesh::MakeCartesian2D(3,
                                     1,
                                     mfem::Element::Type::QUADRILATERAL,
                                     true,
                                     3.0,
                                     1.0));

      for (int i = 0; i < mesh->GetNBE(); i++)
      {
         mfem::Element *be = mesh->GetBdrElement(i);
         mfem::Array<int> vertices;
         be->GetVertices(vertices);

         mfem::real_t *coords1 = mesh->GetVertex(vertices[0]);
         mfem::real_t *coords2 = mesh->GetVertex(vertices[1]);
         mfem::Vector center(2);
         center(0) = 0.5 * (coords1[0] + coords2[0]);
         center(1) = 0.5 * (coords1[1] + coords2[1]);

         if (std::abs(center(0) - 0.0) < 1.0e-10)
         {
            be->SetAttribute(1);
         }
         else
         {
            be->SetAttribute(2);
         }
      }
      mesh->SetAttributes();
      return reinterpret_cast<CMFEM_Mesh *>(mesh);
   }

   void *CMFEM_Ex37FilterSolver_NewMeshOrderEpsAi(CMFEM_Mesh *mesh,
                                                  int order,
                                                  double epsilon_sq,
                                                  const CMFEM_ArrayInt *ess_bdr)
   {
      const mfem::Array<int> *ess_ptr = ess_bdr == nullptr ?
                                        nullptr :
                                        &cmfem::ArrayIntRef(ess_bdr);
      return new Ex37FilterSolverContext(
                cmfem::As<mfem::Mesh>(mesh),
                order,
                static_cast<mfem::real_t>(epsilon_sq),
                ess_ptr);
   }

   void CMFEM_Ex37FilterSolver_Delete(void *solver)
   {
      delete static_cast<Ex37FilterSolverContext *>(solver);
   }

   void CMFEM_Ex37FilterSolver_SolveSigmoidGf(void *solver,
                                              const CMFEM_GridFunction *psi,
                                              CMFEM_GridFunction *rho_filter)
   {
      static_cast<Ex37FilterSolverContext *>(solver)->SolveSigmoid(
         *cmfem::As<const mfem::GridFunction>(psi),
         *cmfem::As<mfem::GridFunction>(rho_filter));
   }

   void CMFEM_Ex37FilterSolver_SolveSedGfGf(void *solver,
                                            double lambda,
                                            double mu,
                                            const CMFEM_GridFunction *u,
                                            const CMFEM_GridFunction *rho_filter,
                                            double rho_min,
                                            CMFEM_GridFunction *w_filter)
   {
      static_cast<Ex37FilterSolverContext *>(solver)->SolveStrainEnergy(
         static_cast<mfem::real_t>(lambda),
         static_cast<mfem::real_t>(mu),
         *cmfem::As<const mfem::GridFunction>(u),
         *cmfem::As<const mfem::GridFunction>(rho_filter),
         static_cast<mfem::real_t>(rho_min),
         *cmfem::As<mfem::GridFunction>(w_filter));
   }

   void *CMFEM_Ex37ElasticitySolver_NewMeshOrderAiCircle(
      CMFEM_Mesh *mesh,
      int order,
      const CMFEM_ArrayInt *ess_bdr,
      double radius,
      const double *center_xy,
      const double *force_xy)
   {
      const mfem::Array<int> *ess_ptr = ess_bdr == nullptr ?
                                        nullptr :
                                        &cmfem::ArrayIntRef(ess_bdr);
      return new Ex37ElasticitySolverContext(
                cmfem::As<mfem::Mesh>(mesh),
                order,
                ess_ptr,
                static_cast<mfem::real_t>(radius),
                center_xy,
                force_xy);
   }

   void CMFEM_Ex37ElasticitySolver_Delete(void *solver)
   {
      delete static_cast<Ex37ElasticitySolverContext *>(solver);
   }

   void CMFEM_Ex37ElasticitySolver_SolveGf(void *solver,
                                           double lambda,
                                           double mu,
                                           const CMFEM_GridFunction *rho_filter,
                                           double rho_min,
                                           CMFEM_GridFunction *u)
   {
      static_cast<Ex37ElasticitySolverContext *>(solver)->Solve(
         static_cast<mfem::real_t>(lambda),
         static_cast<mfem::real_t>(mu),
         *cmfem::As<const mfem::GridFunction>(rho_filter),
         static_cast<mfem::real_t>(rho_min),
         *cmfem::As<mfem::GridFunction>(u));
   }

   double CMFEM_Ex37ElasticitySolver_ComplianceGf(void *solver,
                                                  const CMFEM_GridFunction *u)
   {
      return static_cast<double>(
                static_cast<Ex37ElasticitySolverContext *>(solver)->Compliance(
                   *cmfem::As<const mfem::GridFunction>(u)));
   }

   double CMFEM_Ex37_ProjectSigmoidVolumeGf(CMFEM_GridFunction *psi,
                                            double target_volume,
                                            double tol,
                                            int max_its)
   {
      return static_cast<double>(
                Ex37ProjectVolume(*cmfem::As<mfem::GridFunction>(psi),
                                  static_cast<mfem::real_t>(target_volume),
                                  static_cast<mfem::real_t>(tol),
                                  max_its));
   }

   double CMFEM_Ex37_SigmoidDiffL1GfGf(const CMFEM_GridFunction *psi,
                                       const CMFEM_GridFunction *psi_old)
   {
      mfem::DiffMappedGridFunctionCoefficient diff(
         cmfem::As<const mfem::GridFunction>(psi),
         cmfem::As<const mfem::GridFunction>(psi_old),
         mfem::sigmoid);
      mfem::GridFunction zero(const_cast<mfem::GridFunction *>(cmfem::As<
                                                               const mfem::GridFunction>(psi))->FESpace());
      zero = 0.0;
      return static_cast<double>(zero.ComputeL1Error(diff));
   }

   void CMFEM_Ex37_ProjectSigmoidGf(const CMFEM_GridFunction *psi,
                                    CMFEM_GridFunction *rho_gf)
   {
      mfem::MappedGridFunctionCoefficient rho(
         cmfem::As<const mfem::GridFunction>(psi), mfem::sigmoid);
      cmfem::As<mfem::GridFunction>(rho_gf)->ProjectCoefficient(rho);
   }

   void CMFEM_Ex37_ProjectSimpGf(const CMFEM_GridFunction *rho_filter,
                                 double rho_min,
                                 CMFEM_GridFunction *rho_gf)
   {
      mfem::SIMPInterpolationCoefficient simp(
         const_cast<mfem::GridFunction *>(cmfem::As<const mfem::GridFunction>(
                                             rho_filter)),
         static_cast<mfem::real_t>(rho_min),
         1.0);
      cmfem::As<mfem::GridFunction>(rho_gf)->ProjectCoefficient(simp);
   }

} // extern "C"
