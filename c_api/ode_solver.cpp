#include "common.hpp"

#include <memory>

namespace
{

struct ODESolverHolder
{
   std::unique_ptr<mfem::ODESolver> solver;
};

CMFEM_ASSERT_TYPE(CMFEM_ODESolver, ODESolverHolder);

} // namespace

extern "C" {

   CMFEM_ODESolver *CMFEM_ODESolver_NewType(int type)
   {
      auto solver = mfem::ODESolver::Select(type);
      if (!solver)
      {
         return nullptr;
      }
      return reinterpret_cast<CMFEM_ODESolver *>(
                new ODESolverHolder
      {
         std::move(solver),
      });
   }

   CMFEM_ODESolver *CMFEM_ODESolver_NewImexType(int type)
   {
      auto solver = mfem::ODESolver::SelectIMEX(type);
      if (!solver)
      {
         return nullptr;
      }
      return reinterpret_cast<CMFEM_ODESolver *>(
                new ODESolverHolder
      {
         std::move(solver),
      });
   }

   void CMFEM_ODESolver_Delete(CMFEM_ODESolver *solver)
   {
      delete cmfem::As<ODESolverHolder>(solver);
   }

   void CMFEM_ODESolver_Init(CMFEM_ODESolver *solver,
                             CMFEM_TimeDependentOperator *oper)
   {
      cmfem::As<ODESolverHolder>(solver)->solver->Init(
         *cmfem::As<mfem::TimeDependentOperator>(oper));
   }

   void CMFEM_ODESolver_Step(CMFEM_ODESolver *solver,
                             CMFEM_Vector *u,
                             double *time,
                             double *dt)
   {
      auto t = static_cast<mfem::real_t>(*time);
      auto dt_value = static_cast<mfem::real_t>(*dt);
      cmfem::As<ODESolverHolder>(solver)->solver->Step(cmfem::VectorRef(u),
                                                       t,
                                                       dt_value);
      *time = t;
      *dt = dt_value;
   }

} // extern "C"
