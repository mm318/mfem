#include "diffusion_multigrid.h"
#include "common.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace
{

class DiffusionMultigridImpl : public mfem::GeometricMultigrid
{
private:
   mfem::ConstantCoefficient coeff;

   void ConstructBilinearForm(mfem::FiniteElementSpace &fespace)
   {
      auto *form = new mfem::BilinearForm(&fespace);
      form->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
      form->AddDomainIntegrator(new mfem::DiffusionIntegrator(coeff));
      form->Assemble();
      bfs.Append(form);
   }

   void ConstructCoarseOperatorAndSolver(mfem::FiniteElementSpace &coarse_fespace)
   {
      ConstructBilinearForm(coarse_fespace);

      mfem::OperatorPtr opr;
      opr.SetType(mfem::Operator::ANY_TYPE);
      bfs[0]->FormSystemMatrix(*essentialTrueDofs[0], opr);
      opr.SetOperatorOwner(false);

      auto *pcg = new mfem::CGSolver();
      pcg->SetPrintLevel(-1);
      pcg->SetMaxIter(200);
      pcg->SetRelTol(std::sqrt(1e-4));
      pcg->SetAbsTol(0.0);
      pcg->SetOperator(*opr.Ptr());

      AddLevel(opr.Ptr(), pcg, true, true);
   }

   void ConstructOperatorAndSmoother(mfem::FiniteElementSpace &fespace,
                                     int level)
   {
      const mfem::Array<int> &ess_tdof_list = *essentialTrueDofs[level];
      ConstructBilinearForm(fespace);

      mfem::OperatorPtr opr;
      opr.SetType(mfem::Operator::ANY_TYPE);
      bfs[level]->FormSystemMatrix(ess_tdof_list, opr);
      opr.SetOperatorOwner(false);

      mfem::Vector diag(fespace.GetTrueVSize());
      bfs[level]->AssembleDiagonal(diag);

      auto *smoother = new mfem::OperatorChebyshevSmoother(*opr,
                                                           diag,
                                                           ess_tdof_list,
                                                           2);
      AddLevel(opr.Ptr(), smoother, true, true);
   }

public:
   DiffusionMultigridImpl(mfem::FiniteElementSpaceHierarchy &fespaces,
                          mfem::Array<int> &ess_bdr)
      : mfem::GeometricMultigrid(fespaces, ess_bdr),
        coeff(1.0)
   {
      ConstructCoarseOperatorAndSolver(fespaces.GetFESpaceAtLevel(0));
      for (int level = 1; level < fespaces.GetNumLevels(); ++level)
      {
         ConstructOperatorAndSmoother(fespaces.GetFESpaceAtLevel(level), level);
      }
   }
};

class DiffusionMultigridContext
{
private:
   std::vector<std::unique_ptr<mfem::FiniteElementCollection> > collections;
   std::unique_ptr<mfem::FiniteElementSpaceHierarchy> fespaces;
   std::unique_ptr<DiffusionMultigridImpl> multigrid;

public:
   DiffusionMultigridContext(mfem::Mesh *mesh,
                             int geometric_refinements,
                             int order_refinements,
                             const mfem::Array<int> &ess_bdr)
   {
      const int dim = mesh->Dimension();
      auto *fec = new mfem::H1_FECollection(1, dim);
      collections.emplace_back(fec);
      auto *coarse_fespace = new mfem::FiniteElementSpace(mesh, fec);
      fespaces = std::make_unique<mfem::FiniteElementSpaceHierarchy>(
                    mesh,
                    coarse_fespace,
                    true,
                    true);

      for (int level = 0; level < geometric_refinements; ++level)
      {
         fespaces->AddUniformlyRefinedLevel();
      }

      for (int level = 0; level < order_refinements; ++level)
      {
         auto *next_fec = new mfem::H1_FECollection(static_cast<int>(
                                                       std::pow(2.0, level + 1.0)),
                                                    dim);
         collections.emplace_back(next_fec);
         fespaces->AddOrderRefinedLevel(next_fec);
      }

      mfem::Array<int> ess_bdr_copy(ess_bdr);
      multigrid = std::make_unique<DiffusionMultigridImpl>(*fespaces,
                                                           ess_bdr_copy);
   }

   mfem::FiniteElementSpace &GetFinestFESpace()
   {
      return fespaces->GetFinestFESpace();
   }

   mfem::Mesh &GetFinestMesh()
   {
      return *fespaces->GetFinestFESpace().GetMesh();
   }

   DiffusionMultigridImpl &Preconditioner()
   {
      return *multigrid;
   }

   void SetVCycle(int pre_smoothing, int post_smoothing)
   {
      multigrid->SetCycleType(mfem::Multigrid::CycleType::VCYCLE,
                              pre_smoothing,
                              post_smoothing);
   }

   void FormFineLinearSystem(mfem::GridFunction &x,
                             mfem::LinearForm &b,
                             mfem::OperatorPtr &A,
                             mfem::Vector &X,
                             mfem::Vector &B)
   {
      multigrid->FormFineLinearSystem(x, b, A, X, B);
   }

   void RecoverFineFEMSolution(const mfem::Vector &X,
                               const mfem::LinearForm &b,
                               mfem::GridFunction &x)
   {
      multigrid->RecoverFineFEMSolution(X, b, x);
   }
};

CMFEM_ASSERT_TYPE(CMFEM_DiffusionMultigrid, DiffusionMultigridContext);

} // namespace

extern "C" {

   CMFEM_DiffusionMultigrid *CMFEM_DiffusionMultigrid_NewMeshGrOrAi(
      CMFEM_Mesh *mesh,
      int geometric_refinements,
      int order_refinements,
      const CMFEM_ArrayInt *ess_bdr)
   {
      return reinterpret_cast<CMFEM_DiffusionMultigrid *>(
                new DiffusionMultigridContext(
                   cmfem::As<mfem::Mesh>(mesh),
                   geometric_refinements,
                   order_refinements,
                   cmfem::ArrayIntRef(ess_bdr)));
   }

   void CMFEM_DiffusionMultigrid_Delete(CMFEM_DiffusionMultigrid *multigrid)
   {
      delete cmfem::As<DiffusionMultigridContext>(multigrid);
   }

   CMFEM_FiniteElementSpace *CMFEM_DiffusionMultigrid_GetFinestFESpace(
      CMFEM_DiffusionMultigrid *multigrid)
   {
      return reinterpret_cast<CMFEM_FiniteElementSpace *>(
                &cmfem::As<DiffusionMultigridContext>(multigrid)->GetFinestFESpace());
   }

   CMFEM_Mesh *CMFEM_DiffusionMultigrid_GetFinestMesh(
      CMFEM_DiffusionMultigrid *multigrid)
   {
      return reinterpret_cast<CMFEM_Mesh *>(
                &cmfem::As<DiffusionMultigridContext>(multigrid)->GetFinestMesh());
   }

   void CMFEM_DiffusionMultigrid_SetCycleTypeV(
      CMFEM_DiffusionMultigrid *multigrid,
      int pre_smoothing,
      int post_smoothing)
   {
      cmfem::As<DiffusionMultigridContext>(multigrid)->SetVCycle(pre_smoothing,
                                                                 post_smoothing);
   }

   void CMFEM_DiffusionMultigrid_FormFineLinearSystemOp(
      CMFEM_DiffusionMultigrid *multigrid,
      CMFEM_GridFunction *x,
      CMFEM_LinearForm *b,
      CMFEM_OperatorPtr *A,
      CMFEM_Vector *X,
      CMFEM_Vector *B)
   {
      cmfem::As<DiffusionMultigridContext>(multigrid)->FormFineLinearSystem(
         *cmfem::As<mfem::GridFunction>(x),
         *cmfem::As<mfem::LinearForm>(b),
         cmfem::OperatorPtrRef(A),
         cmfem::VectorRef(X),
         cmfem::VectorRef(B));
   }

   void CMFEM_DiffusionMultigrid_RecoverFineFEMSolution(
      CMFEM_DiffusionMultigrid *multigrid,
      const CMFEM_Vector *X,
      const CMFEM_LinearForm *b,
      CMFEM_GridFunction *x)
   {
      cmfem::As<DiffusionMultigridContext>(multigrid)->RecoverFineFEMSolution(
         cmfem::VectorRef(X),
         *cmfem::As<const mfem::LinearForm>(b),
         *cmfem::As<mfem::GridFunction>(x));
   }

   void CMFEM_DiffusionMultigrid_PCGOp(
      CMFEM_DiffusionMultigrid *multigrid,
      const CMFEM_OperatorPtr *A,
      const CMFEM_Vector *B,
      CMFEM_Vector *X,
      int print_iter,
      int max_iter,
      double rtol,
      double atol)
   {
      mfem::PCG(*cmfem::OperatorPtrRef(A),
                cmfem::As<DiffusionMultigridContext>(multigrid)->Preconditioner(),
                cmfem::VectorRef(B),
                cmfem::VectorRef(X),
                print_iter,
                max_iter,
                static_cast<mfem::real_t>(rtol),
                static_cast<mfem::real_t>(atol));
   }

} // extern "C"
