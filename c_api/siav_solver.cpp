#include "common.hpp"

#include <memory>
#include <unordered_map>

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_SIAVSolver, mfem::SIAVSolver);

class COperatorAdapter : public mfem::Operator
{
private:
   CMFEM_OperatorMultCallback callback;
   void *context;

public:
   COperatorAdapter(int size,
                    CMFEM_OperatorMultCallback callback_,
                    void *context_)
      : mfem::Operator(size),
        callback(callback_),
        context(context_) { }

   void Mult(const mfem::Vector &x, mfem::Vector &y) const override
   {
      callback(reinterpret_cast<const CMFEM_Vector *>(&x),
               reinterpret_cast<CMFEM_Vector *>(&y),
               context);
   }
};

class CTimeDependentOperatorAdapter : public mfem::TimeDependentOperator
{
private:
   CMFEM_TimeDependentOperatorMultCallback callback;
   void *context;

public:
   CTimeDependentOperatorAdapter(int size,
                                 CMFEM_TimeDependentOperatorMultCallback callback_,
                                 void *context_)
      : mfem::TimeDependentOperator(size),
        callback(callback_),
        context(context_) { }

   void Mult(const mfem::Vector &x, mfem::Vector &y) const override
   {
      callback(reinterpret_cast<const CMFEM_Vector *>(&x),
               reinterpret_cast<CMFEM_Vector *>(&y),
               GetTime(),
               context);
   }
};

struct SIAVContext
{
   std::unique_ptr<COperatorAdapter> p;
   std::unique_ptr<CTimeDependentOperatorAdapter> f;
};

std::unordered_map<const CMFEM_SIAVSolver *, SIAVContext> &Contexts()
{
   static std::unordered_map<const CMFEM_SIAVSolver *, SIAVContext> contexts;
   return contexts;
}

} // namespace

extern "C" {

   CMFEM_SIAVSolver *CMFEM_SIAVSolver_NewOrder(int order)
   {
      return reinterpret_cast<CMFEM_SIAVSolver *>(new mfem::SIAVSolver(order));
   }

   void CMFEM_SIAVSolver_Delete(CMFEM_SIAVSolver *solver)
   {
      Contexts().erase(solver);
      delete cmfem::As<mfem::SIAVSolver>(solver);
   }

   void CMFEM_SIAVSolver_Init(CMFEM_SIAVSolver *solver,
                              int size,
                              CMFEM_OperatorMultCallback p_mult,
                              void *p_context,
                              CMFEM_TimeDependentOperatorMultCallback f_mult,
                              void *f_context)
   {
      auto context = SIAVContext{
         std::make_unique<COperatorAdapter>(size, p_mult, p_context),
         std::make_unique<CTimeDependentOperatorAdapter>(size, f_mult, f_context),
      };

      auto &contexts = Contexts();
      contexts[solver] = std::move(context);
      cmfem::As<mfem::SIAVSolver>(solver)->Init(*contexts[solver].p,
                                                *contexts[solver].f);
   }

   void CMFEM_SIAVSolver_Step(CMFEM_SIAVSolver *solver,
                              CMFEM_Vector *q,
                              CMFEM_Vector *p,
                              double *time,
                              double *dt)
   {
      auto t = static_cast<mfem::real_t>(*time);
      auto dt_value = static_cast<mfem::real_t>(*dt);
      cmfem::As<mfem::SIAVSolver>(solver)->Step(cmfem::VectorRef(q),
                                                cmfem::VectorRef(p),
                                                t,
                                                dt_value);
      *time = t;
      *dt = dt_value;
   }

} // extern "C"
