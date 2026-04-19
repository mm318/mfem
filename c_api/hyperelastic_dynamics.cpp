#include "common.hpp"

namespace
{

class Ex10ReducedSystemOperator;

class Ex10HyperelasticOperator
{
private:
   mfem::FiniteElementSpace &fespace;
   mfem::BilinearForm M;
   mfem::BilinearForm S;
   mfem::NonlinearForm H;
   mfem::real_t viscosity;
   mfem::HyperelasticModel *model;
   mfem::CGSolver M_solver;
   mfem::DSmoother M_prec;
   Ex10ReducedSystemOperator *reduced_oper;
   mfem::NewtonSolver newton_solver;
   mfem::Solver *J_solver;
   mfem::Solver *J_prec;
   mutable mfem::Vector z;

public:
   Ex10HyperelasticOperator(mfem::FiniteElementSpace &f,
                            const mfem::Array<int> &ess_bdr,
                            mfem::real_t visc,
                            mfem::real_t mu,
                            mfem::real_t K);

   int StateSize() const { return 2 * fespace.GetTrueVSize(); }

   void SetState(const mfem::GridFunction &velocity,
                 const mfem::GridFunction &position,
                 mfem::Vector &state) const;
   void SetFields(const mfem::Vector &state,
                  mfem::GridFunction &velocity,
                  mfem::GridFunction &position) const;
   void Mult(const mfem::Vector &state, mfem::Vector &dstate_dt) const;
   void ImplicitSolve(mfem::real_t dt,
                      const mfem::Vector &state,
                      mfem::Vector &k);
   mfem::real_t ElasticEnergy(const mfem::Vector &state) const;
   mfem::real_t KineticEnergy(const mfem::Vector &state) const;
   void GetElasticEnergyDensity(const mfem::Vector &state,
                                mfem::GridFunction &w) const;

   ~Ex10HyperelasticOperator();
};

class Ex10ReducedSystemOperator : public mfem::Operator
{
private:
   mfem::BilinearForm *M;
   mfem::BilinearForm *S;
   mfem::NonlinearForm *H;
   mutable mfem::SparseMatrix *Jacobian;
   mfem::real_t dt;
   const mfem::Vector *v;
   const mfem::Vector *x;
   mutable mfem::Vector w;
   mutable mfem::Vector z;

public:
   Ex10ReducedSystemOperator(mfem::BilinearForm *M_,
                             mfem::BilinearForm *S_,
                             mfem::NonlinearForm *H_)
      : mfem::Operator(M_->Height()),
        M(M_),
        S(S_),
        H(H_),
        Jacobian(nullptr),
        dt(0.0),
        v(nullptr),
        x(nullptr),
        w(height),
        z(height) { }

   void SetParameters(mfem::real_t dt_,
                      const mfem::Vector *v_,
                      const mfem::Vector *x_)
   {
      dt = dt_;
      v = v_;
      x = x_;
   }

   void Mult(const mfem::Vector &k, mfem::Vector &y) const override
   {
      add(*v, dt, k, w);
      add(*x, dt, w, z);
      H->Mult(z, y);
      M->AddMult(k, y);
      S->AddMult(w, y);
   }

   mfem::Operator &GetGradient(const mfem::Vector &k) const override
   {
      delete Jacobian;
      Jacobian = Add(1.0, M->SpMat(), dt, S->SpMat());
      add(*v, dt, k, w);
      add(*x, dt, w, z);
      auto *grad_H = dynamic_cast<mfem::SparseMatrix *>(&H->GetGradient(z));
      Jacobian->Add(dt * dt, *grad_H);
      return *Jacobian;
   }

   ~Ex10ReducedSystemOperator() override
   {
      delete Jacobian;
   }
};

class ElasticEnergyCoefficient : public mfem::Coefficient
{
private:
   mfem::HyperelasticModel &model;
   const mfem::GridFunction &x;
   mutable mfem::DenseMatrix J;

public:
   ElasticEnergyCoefficient(mfem::HyperelasticModel &m,
                            const mfem::GridFunction &x_)
      : model(m), x(x_) { }

   mfem::real_t Eval(mfem::ElementTransformation &T,
                     const mfem::IntegrationPoint &ip) override
   {
      model.SetTransformation(T);
      x.GetVectorGradient(T, J);
      return model.EvalW(J) / J.Det();
   }
};

void SplitState(const mfem::Vector &state, mfem::Vector &velocity,
                mfem::Vector &position)
{
   const int size = state.Size() / 2;
   velocity.NewDataAndSize(const_cast<mfem::real_t *>(state.GetData()), size);
   position.NewDataAndSize(const_cast<mfem::real_t *>(state.GetData()) + size,
                           size);
}

Ex10HyperelasticOperator::Ex10HyperelasticOperator(
   mfem::FiniteElementSpace &f,
   const mfem::Array<int> &ess_bdr,
   mfem::real_t visc,
   mfem::real_t mu,
   mfem::real_t K)
   : fespace(f),
     M(&fespace),
     S(&fespace),
     H(&fespace),
     viscosity(visc),
     model(nullptr),
     reduced_oper(nullptr),
     J_solver(nullptr),
     J_prec(nullptr),
     z(f.GetTrueVSize())
{
#if defined(MFEM_USE_DOUBLE)
   const mfem::real_t rel_tol = 1e-8;
   const mfem::real_t newton_abs_tol = 0.0;
#elif defined(MFEM_USE_SINGLE)
   const mfem::real_t rel_tol = 1e-3;
   const mfem::real_t newton_abs_tol = 1e-4;
#else
#error "Only single and double precision are supported!"
#endif
   const int skip_zero_entries = 0;

   const mfem::real_t ref_density = 1.0;
   mfem::ConstantCoefficient rho0(ref_density);
   M.AddDomainIntegrator(new mfem::VectorMassIntegrator(rho0));
   M.Assemble(skip_zero_entries);

   mfem::Array<int> ess_tdof_list;
   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   mfem::SparseMatrix tmp;
   M.FormSystemMatrix(ess_tdof_list, tmp);

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(30);
   M_solver.SetPrintLevel(0);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(M.SpMat());

   model = new mfem::NeoHookeanModel(mu, K);
   H.AddDomainIntegrator(new mfem::HyperelasticNLFIntegrator(model));
   H.SetEssentialTrueDofs(ess_tdof_list);

   mfem::ConstantCoefficient visc_coeff(viscosity);
   S.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(visc_coeff));
   S.Assemble(skip_zero_entries);
   S.FormSystemMatrix(ess_tdof_list, tmp);

   reduced_oper = new Ex10ReducedSystemOperator(&M, &S, &H);

#ifndef MFEM_USE_SUITESPARSE
   J_prec = new mfem::DSmoother(1);
   auto *J_minres = new mfem::MINRESSolver;
   J_minres->SetRelTol(rel_tol);
   J_minres->SetAbsTol(0.0);
   J_minres->SetMaxIter(300);
   J_minres->SetPrintLevel(-1);
   J_minres->SetPreconditioner(*J_prec);
   J_solver = J_minres;
#else
   J_solver = new mfem::UMFPackSolver;
   J_prec = nullptr;
#endif

   newton_solver.iterative_mode = false;
   newton_solver.SetSolver(*J_solver);
   newton_solver.SetOperator(*reduced_oper);
   newton_solver.SetPrintLevel(1);
   newton_solver.SetRelTol(rel_tol);
   newton_solver.SetAbsTol(newton_abs_tol);
   newton_solver.SetMaxIter(10);
}

void Ex10HyperelasticOperator::SetState(const mfem::GridFunction &velocity,
                                        const mfem::GridFunction &position,
                                        mfem::Vector &state) const
{
   state.SetSize(StateSize());
   mfem::Vector velocity_state(state.GetData(), z.Size());
   mfem::Vector position_state(state.GetData() + z.Size(), z.Size());
   velocity.GetTrueDofs(velocity_state);
   position.GetTrueDofs(position_state);
}

void Ex10HyperelasticOperator::SetFields(const mfem::Vector &state,
                                         mfem::GridFunction &velocity,
                                         mfem::GridFunction &position) const
{
   mfem::Vector velocity_state;
   mfem::Vector position_state;
   SplitState(state, velocity_state, position_state);
   velocity.SetFromTrueDofs(velocity_state);
   position.SetFromTrueDofs(position_state);
}

void Ex10HyperelasticOperator::Mult(const mfem::Vector &state,
                                    mfem::Vector &dstate_dt) const
{
   dstate_dt.SetSize(StateSize());
   mfem::Vector velocity;
   mfem::Vector position;
   mfem::Vector dvelocity_dt(dstate_dt.GetData(), z.Size());
   mfem::Vector dposition_dt(dstate_dt.GetData() + z.Size(), z.Size());
   SplitState(state, velocity, position);

   H.Mult(position, z);
   if (viscosity != 0.0)
   {
      S.AddMult(velocity, z);
   }
   z.Neg();
   M_solver.Mult(z, dvelocity_dt);
   dposition_dt = velocity;
}

void Ex10HyperelasticOperator::ImplicitSolve(mfem::real_t dt,
                                             const mfem::Vector &state,
                                             mfem::Vector &k)
{
   k.SetSize(StateSize());
   mfem::Vector velocity;
   mfem::Vector position;
   mfem::Vector dvelocity_dt(k.GetData(), z.Size());
   mfem::Vector dposition_dt(k.GetData() + z.Size(), z.Size());
   SplitState(state, velocity, position);

   reduced_oper->SetParameters(dt, &velocity, &position);
   mfem::Vector zero;
   newton_solver.Mult(zero, dvelocity_dt);
   MFEM_VERIFY(newton_solver.GetConverged(), "Newton solver did not converge.");
   add(velocity, dt, dvelocity_dt, dposition_dt);
}

mfem::real_t Ex10HyperelasticOperator::ElasticEnergy(const mfem::Vector &state)
const
{
   mfem::Vector velocity;
   mfem::Vector position;
   SplitState(state, velocity, position);
   return H.GetEnergy(position);
}

mfem::real_t Ex10HyperelasticOperator::KineticEnergy(const mfem::Vector &state)
const
{
   mfem::Vector velocity;
   mfem::Vector position;
   SplitState(state, velocity, position);
   return 0.5 * M.InnerProduct(velocity, velocity);
}

void Ex10HyperelasticOperator::GetElasticEnergyDensity(const mfem::Vector
                                                       &state,
                                                       mfem::GridFunction &w) const
{
   mfem::Vector velocity;
   mfem::Vector position_state;
   SplitState(state, velocity, position_state);
   mfem::GridFunction position(&fespace);
   position.SetFromTrueDofs(position_state);
   ElasticEnergyCoefficient w_coeff(*model, position);
   w.ProjectCoefficient(w_coeff);
}

Ex10HyperelasticOperator::~Ex10HyperelasticOperator()
{
   delete J_solver;
   delete J_prec;
   delete reduced_oper;
   delete model;
}

} // namespace

extern "C" {

   void *CMFEM_Ex10HyperelasticOperator_New(CMFEM_FiniteElementSpace *fespace,
                                            const CMFEM_ArrayInt *ess_bdr,
                                            double viscosity,
                                            double mu,
                                            double bulk_modulus)
   {
      return new Ex10HyperelasticOperator(
                *cmfem::As<mfem::FiniteElementSpace>(fespace),
                cmfem::ArrayIntRef(const_cast<CMFEM_ArrayInt *>(ess_bdr)),
                static_cast<mfem::real_t>(viscosity),
                static_cast<mfem::real_t>(mu),
                static_cast<mfem::real_t>(bulk_modulus));
   }

   void CMFEM_Ex10HyperelasticOperator_Delete(void *oper)
   {
      delete static_cast<Ex10HyperelasticOperator *>(oper);
   }

   int CMFEM_Ex10HyperelasticOperator_StateSize(const void *oper)
   {
      return static_cast<const Ex10HyperelasticOperator *>(oper)->StateSize();
   }

   void CMFEM_Ex10HyperelasticOperator_SetState(void *oper,
                                                const CMFEM_GridFunction *velocity,
                                                const CMFEM_GridFunction *position,
                                                CMFEM_Vector *state)
   {
      static_cast<Ex10HyperelasticOperator *>(oper)->SetState(
         *cmfem::As<const mfem::GridFunction>(velocity),
         *cmfem::As<const mfem::GridFunction>(position),
         cmfem::VectorRef(state));
   }

   void CMFEM_Ex10HyperelasticOperator_SetFields(const void *oper,
                                                 const CMFEM_Vector *state,
                                                 CMFEM_GridFunction *velocity,
                                                 CMFEM_GridFunction *position)
   {
      static_cast<const Ex10HyperelasticOperator *>(oper)->SetFields(
         cmfem::VectorRef(const_cast<CMFEM_Vector *>(state)),
         *cmfem::As<mfem::GridFunction>(velocity),
         *cmfem::As<mfem::GridFunction>(position));
   }

   void CMFEM_Ex10HyperelasticOperator_Mult(void *oper,
                                            const CMFEM_Vector *state,
                                            CMFEM_Vector *dstate_dt)
   {
      static_cast<Ex10HyperelasticOperator *>(oper)->Mult(
         cmfem::VectorRef(const_cast<CMFEM_Vector *>(state)),
         cmfem::VectorRef(dstate_dt));
   }

   void CMFEM_Ex10HyperelasticOperator_ImplicitSolve(void *oper,
                                                     double dt,
                                                     const CMFEM_Vector *state,
                                                     CMFEM_Vector *k)
   {
      static_cast<Ex10HyperelasticOperator *>(oper)->ImplicitSolve(
         static_cast<mfem::real_t>(dt),
         cmfem::VectorRef(const_cast<CMFEM_Vector *>(state)),
         cmfem::VectorRef(k));
   }

   double CMFEM_Ex10HyperelasticOperator_ElasticEnergy(void *oper,
                                                       const CMFEM_Vector *state)
   {
      return static_cast<double>(
                static_cast<Ex10HyperelasticOperator *>(oper)->ElasticEnergy(
                   cmfem::VectorRef(const_cast<CMFEM_Vector *>(state))));
   }

   double CMFEM_Ex10HyperelasticOperator_KineticEnergy(void *oper,
                                                       const CMFEM_Vector *state)
   {
      return static_cast<double>(
                static_cast<Ex10HyperelasticOperator *>(oper)->KineticEnergy(
                   cmfem::VectorRef(const_cast<CMFEM_Vector *>(state))));
   }

   void CMFEM_Ex10HyperelasticOperator_GetElasticEnergyDensity(
      void *oper,
      const CMFEM_Vector *state,
      CMFEM_GridFunction *energy_density)
   {
      static_cast<Ex10HyperelasticOperator *>(oper)->GetElasticEnergyDensity(
         cmfem::VectorRef(const_cast<CMFEM_Vector *>(state)),
         *cmfem::As<mfem::GridFunction>(energy_density));
   }

} // extern "C"
