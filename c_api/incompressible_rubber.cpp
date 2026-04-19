#include "common.hpp"

#include <iomanip>
#include <memory>
#include <string>

namespace
{

class GeneralResidualMonitor : public mfem::IterativeSolverMonitor
{
private:
   const std::string prefix;
   int print_level;
   mutable mfem::real_t norm0;

public:
   GeneralResidualMonitor(const std::string &prefix_, int print_level_)
      : prefix(prefix_), print_level(print_level_), norm0(0.0) { }

   void MonitorResidual(int it,
                        mfem::real_t norm,
                        const mfem::Vector &r,
                        bool final) override
   {
      (void)r;
      if ((print_level == 1 && !final) ||
          (print_level == 3 && (final || it == 0)))
      {
         mfem::out << prefix << " iteration " << std::setw(2) << it
                   << " : ||r|| = " << norm;
         if (it > 0)
         {
            mfem::out << ",  ||r||/||r_0|| = " << norm / norm0;
         }
         else
         {
            norm0 = norm;
         }
         mfem::out << '\n';
      }
   }
};

class JacobianPreconditioner : public mfem::Solver
{
private:
   mfem::Array<mfem::FiniteElementSpace *> spaces;
   mfem::Array<int> &block_true_offsets;
   mfem::BlockOperator *jacobian;
   mfem::real_t gamma;
   mfem::SparseMatrix *pressure_mass;
   mfem::Solver *mass_pcg;
   mfem::Solver *mass_prec;
   mfem::Solver *stiff_pcg;
   mfem::Solver *stiff_prec;

public:
   JacobianPreconditioner(mfem::Array<mfem::FiniteElementSpace *> &fes,
                          mfem::SparseMatrix &mass,
                          mfem::Array<int> &offsets)
      : mfem::Solver(offsets[2]),
        block_true_offsets(offsets),
        jacobian(nullptr),
        gamma(0.00001),
        pressure_mass(&mass),
        mass_pcg(nullptr),
        mass_prec(nullptr),
        stiff_pcg(nullptr),
        stiff_prec(nullptr)
   {
      fes.Copy(spaces);

      auto *mass_prec_gs = new mfem::GSSmoother(*pressure_mass);
      mass_prec = mass_prec_gs;

      auto *mass_pcg_iter = new mfem::CGSolver();
      mass_pcg_iter->SetRelTol(1e-12);
      mass_pcg_iter->SetAbsTol(1e-12);
      mass_pcg_iter->SetMaxIter(200);
      mass_pcg_iter->SetPrintLevel(0);
      mass_pcg_iter->SetPreconditioner(*mass_prec);
      mass_pcg_iter->SetOperator(*pressure_mass);
      mass_pcg_iter->iterative_mode = false;
      mass_pcg = mass_pcg_iter;
   }

   void Mult(const mfem::Vector &k, mfem::Vector &y) const override
   {
      mfem::Vector disp_in(k.GetData() + block_true_offsets[0],
                           block_true_offsets[1] - block_true_offsets[0]);
      mfem::Vector pres_in(k.GetData() + block_true_offsets[1],
                           block_true_offsets[2] - block_true_offsets[1]);
      mfem::Vector disp_out(y.GetData() + block_true_offsets[0],
                            block_true_offsets[1] - block_true_offsets[0]);
      mfem::Vector pres_out(y.GetData() + block_true_offsets[1],
                            block_true_offsets[2] - block_true_offsets[1]);
      mfem::Vector temp(block_true_offsets[1] - block_true_offsets[0]);
      mfem::Vector temp2(block_true_offsets[1] - block_true_offsets[0]);

      mass_pcg->Mult(pres_in, pres_out);
      pres_out *= -gamma;

      jacobian->GetBlock(0, 1).Mult(pres_out, temp);
      subtract(disp_in, temp, temp2);

      stiff_pcg->Mult(temp2, disp_out);
   }

   void SetOperator(const mfem::Operator &op) override
   {
      jacobian = (mfem::BlockOperator *)&op;
      if (stiff_prec == nullptr)
      {
         auto *stiff_prec_gs = new mfem::GSSmoother();
         stiff_prec = stiff_prec_gs;

         auto *stiff_pcg_iter = new mfem::GMRESSolver();
         stiff_pcg_iter->SetRelTol(1e-8);
         stiff_pcg_iter->SetAbsTol(1e-8);
         stiff_pcg_iter->SetMaxIter(200);
         stiff_pcg_iter->SetPrintLevel(0);
         stiff_pcg_iter->SetPreconditioner(*stiff_prec);
         stiff_pcg_iter->iterative_mode = false;
         stiff_pcg = stiff_pcg_iter;
      }

      stiff_pcg->SetOperator(jacobian->GetBlock(0, 0));
   }

   ~JacobianPreconditioner() override
   {
      delete mass_pcg;
      delete mass_prec;
      delete stiff_prec;
      delete stiff_pcg;
   }
};

class RubberOperator : public mfem::Operator
{
private:
   mfem::Array<mfem::FiniteElementSpace *> spaces;
   mfem::BlockNonlinearForm *Hform;
   mfem::SparseMatrix *pressure_mass;
   mutable mfem::NewtonSolver newton_solver;
   mutable GeneralResidualMonitor newton_monitor;
   mfem::Solver *j_solver;
   GeneralResidualMonitor j_monitor;
   mfem::Solver *j_prec;
   mfem::Coefficient &mu;
   mfem::Array<int> &block_true_offsets;

public:
   RubberOperator(mfem::Array<mfem::FiniteElementSpace *> &fes,
                  mfem::Array<mfem::Array<int> *> &ess_bdr,
                  mfem::Array<int> &offsets,
                  mfem::real_t rel_tol,
                  mfem::real_t abs_tol,
                  int iter,
                  mfem::Coefficient &c_mu)
      : mfem::Operator(fes[0]->GetTrueVSize() + fes[1]->GetTrueVSize()),
        Hform(nullptr),
        pressure_mass(nullptr),
        newton_solver(),
        newton_monitor("Newton", 1),
        j_solver(nullptr),
        j_monitor("  GMRES", 3),
        j_prec(nullptr),
        mu(c_mu),
        block_true_offsets(offsets)
   {
      mfem::Array<mfem::Vector *> rhs(2);
      rhs = nullptr;
      fes.Copy(spaces);

      Hform = new mfem::BlockNonlinearForm(spaces);
      Hform->AddDomainIntegrator(new mfem::IncompressibleNeoHookeanIntegrator(mu));
      Hform->SetEssentialBC(ess_bdr, rhs);

      auto *a = new mfem::BilinearForm(spaces[1]);
      mfem::ConstantCoefficient one(1.0);
      a->AddDomainIntegrator(new mfem::MassIntegrator(one));
      a->Assemble();
      a->Finalize();

      mfem::OperatorPtr op;
      mfem::Array<int> p_ess_tdofs;
      a->FormSystemMatrix(p_ess_tdofs, op);
      pressure_mass = a->LoseMat();
      delete a;

      j_prec = new JacobianPreconditioner(fes, *pressure_mass, block_true_offsets);

      auto *j_gmres = new mfem::GMRESSolver();
      j_gmres->iterative_mode = false;
      j_gmres->SetRelTol(1e-12);
      j_gmres->SetAbsTol(1e-12);
      j_gmres->SetMaxIter(300);
      j_gmres->SetPrintLevel(-1);
      j_gmres->SetMonitor(j_monitor);
      j_gmres->SetPreconditioner(*j_prec);
      j_solver = j_gmres;

      newton_solver.iterative_mode = true;
      newton_solver.SetSolver(*j_solver);
      newton_solver.SetOperator(*this);
      newton_solver.SetPrintLevel(-1);
      newton_solver.SetMonitor(newton_monitor);
      newton_solver.SetRelTol(rel_tol);
      newton_solver.SetAbsTol(abs_tol);
      newton_solver.SetMaxIter(iter);
   }

   void Solve(mfem::Vector &xp) const
   {
      mfem::Vector zero;
      newton_solver.Mult(zero, xp);
      MFEM_VERIFY(newton_solver.GetConverged(), "Newton Solver did not converge.");
   }

   void Mult(const mfem::Vector &x, mfem::Vector &y) const override
   {
      Hform->Mult(x, y);
   }

   mfem::Operator &GetGradient(const mfem::Vector &xp) const override
   {
      return Hform->GetGradient(xp);
   }

   ~RubberOperator() override
   {
      delete Hform;
      delete pressure_mass;
      delete j_solver;
      delete j_prec;
   }
};

class Ex19RubberOperatorContext
{
private:
   mfem::Array<mfem::FiniteElementSpace *> spaces;
   mfem::Array<mfem::Array<int> *> ess_bdr;
   mfem::Array<int> ess_bdr_u;
   mfem::Array<int> ess_bdr_p;
   mfem::Array<int> block_true_offsets;
   mfem::ConstantCoefficient mu;
   std::unique_ptr<RubberOperator> oper;

public:
   Ex19RubberOperatorContext(mfem::FiniteElementSpace &displacement_space,
                             mfem::FiniteElementSpace &pressure_space,
                             const mfem::Array<int> &ess_bdr_u_in,
                             const mfem::Array<int> &ess_bdr_p_in,
                             mfem::real_t rel_tol,
                             mfem::real_t abs_tol,
                             int max_iter,
                             mfem::real_t mu_value)
      : spaces(2),
        ess_bdr(2),
        ess_bdr_u(ess_bdr_u_in),
        ess_bdr_p(ess_bdr_p_in),
        block_true_offsets(3),
        mu(mu_value),
        oper()
   {
      spaces[0] = &displacement_space;
      spaces[1] = &pressure_space;
      ess_bdr[0] = &ess_bdr_u;
      ess_bdr[1] = &ess_bdr_p;
      block_true_offsets[0] = 0;
      block_true_offsets[1] = displacement_space.GetTrueVSize();
      block_true_offsets[2] = pressure_space.GetTrueVSize();
      block_true_offsets.PartialSum();
      oper = std::make_unique<RubberOperator>(spaces,
                                              ess_bdr,
                                              block_true_offsets,
                                              rel_tol,
                                              abs_tol,
                                              max_iter,
                                              mu);
   }

   void Solve(mfem::GridFunction &position, mfem::GridFunction &pressure) const
   {
      mfem::Vector state(block_true_offsets.Last());
      mfem::Vector disp(state.GetData() + block_true_offsets[0],
                        block_true_offsets[1] - block_true_offsets[0]);
      mfem::Vector pres(state.GetData() + block_true_offsets[1],
                        block_true_offsets[2] - block_true_offsets[1]);

      position.GetTrueDofs(disp);
      pressure.GetTrueDofs(pres);
      oper->Solve(state);
      position.SetFromTrueDofs(disp);
      pressure.SetFromTrueDofs(pres);
   }
};

} // namespace

extern "C" {

   void *CMFEM_Ex19RubberOperator_New(
      CMFEM_FiniteElementSpace *displacement_space,
      CMFEM_FiniteElementSpace *pressure_space,
      const CMFEM_ArrayInt *ess_bdr_u,
      const CMFEM_ArrayInt *ess_bdr_p,
      double rel_tol,
      double abs_tol,
      int max_iter,
      double mu)
   {
      return new Ex19RubberOperatorContext(
                *cmfem::As<mfem::FiniteElementSpace>(displacement_space),
                *cmfem::As<mfem::FiniteElementSpace>(pressure_space),
                cmfem::ArrayIntRef(const_cast<CMFEM_ArrayInt *>(ess_bdr_u)),
                cmfem::ArrayIntRef(const_cast<CMFEM_ArrayInt *>(ess_bdr_p)),
                static_cast<mfem::real_t>(rel_tol),
                static_cast<mfem::real_t>(abs_tol),
                max_iter,
                static_cast<mfem::real_t>(mu));
   }

   void CMFEM_Ex19RubberOperator_Delete(void *oper)
   {
      delete static_cast<Ex19RubberOperatorContext *>(oper);
   }

   void CMFEM_Ex19RubberOperator_Solve(void *oper,
                                       CMFEM_GridFunction *position,
                                       CMFEM_GridFunction *pressure)
   {
      static_cast<Ex19RubberOperatorContext *>(oper)->Solve(
         *cmfem::As<mfem::GridFunction>(position),
         *cmfem::As<mfem::GridFunction>(pressure));
   }

} // extern "C"
