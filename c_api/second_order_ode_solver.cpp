#include "common.hpp"

#include <memory>

namespace
{

struct SecondOrderODESolverHolder
{
   std::unique_ptr<mfem::SecondOrderODESolver> solver;
};

CMFEM_ASSERT_TYPE(CMFEM_SecondOrderODESolver, SecondOrderODESolverHolder);

} // namespace

extern "C" {

   CMFEM_SecondOrderODESolver *CMFEM_SecondOrderODESolver_NewType(int type)
   {
      auto *solver = mfem::SecondOrderODESolver::Select(type);
      if (!solver)
      {
         return nullptr;
      }
      return reinterpret_cast<CMFEM_SecondOrderODESolver *>(
                new SecondOrderODESolverHolder
      {
         std::unique_ptr<mfem::SecondOrderODESolver>(solver),
      });
   }

   void CMFEM_SecondOrderODESolver_Delete(CMFEM_SecondOrderODESolver *solver)
   {
      delete cmfem::As<SecondOrderODESolverHolder>(solver);
   }

   void CMFEM_SecondOrderODESolver_Init(
      CMFEM_SecondOrderODESolver *solver,
      CMFEM_SecondOrderTimeDependentOperator *oper)
   {
      cmfem::As<SecondOrderODESolverHolder>(solver)->solver->Init(
         *cmfem::As<mfem::SecondOrderTimeDependentOperator>(oper));
   }

   void CMFEM_SecondOrderODESolver_Step(CMFEM_SecondOrderODESolver *solver,
                                        CMFEM_Vector *u,
                                        CMFEM_Vector *du_dt,
                                        double *time,
                                        double *dt)
   {
      auto t = static_cast<mfem::real_t>(*time);
      auto dt_value = static_cast<mfem::real_t>(*dt);
      cmfem::As<SecondOrderODESolverHolder>(solver)->solver->Step(
         cmfem::VectorRef(u),
         cmfem::VectorRef(du_dt),
         t,
         dt_value);
      *time = t;
      *dt = dt_value;
   }

} // extern "C"
