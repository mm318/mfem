#include "common.hpp"

namespace
{

class CTimeDependentOperatorAdapter final : public mfem::TimeDependentOperator
{
private:
   CMFEM_TimeDependentOperatorMultCallback mult;
   void *mult_context;
   CMFEM_TimeDependentOperatorImplicitSolveCallback implicit_solve;
   void *implicit_context;

public:
   CTimeDependentOperatorAdapter(
      int size,
      mfem::TimeDependentOperator::Type type,
      CMFEM_TimeDependentOperatorMultCallback mult_,
      void *mult_context_,
      CMFEM_TimeDependentOperatorImplicitSolveCallback implicit_solve_,
      void *implicit_context_)
      : mfem::TimeDependentOperator(size, 0.0, type),
        mult(mult_),
        mult_context(mult_context_),
        implicit_solve(implicit_solve_),
        implicit_context(implicit_context_) { }

   void Mult(const mfem::Vector &u, mfem::Vector &du_dt) const override
   {
      if (!mult)
      {
         MFEM_ABORT("CMFEM time-dependent operator missing Mult callback");
      }
      mult(reinterpret_cast<const CMFEM_Vector *>(&u),
           reinterpret_cast<CMFEM_Vector *>(&du_dt),
           GetTime(),
           mult_context);
   }

   void ImplicitSolve(const mfem::real_t gamma,
                      const mfem::Vector &u,
                      mfem::Vector &k) override
   {
      if (!implicit_solve)
      {
         MFEM_ABORT(
            "CMFEM time-dependent operator missing ImplicitSolve callback");
      }
      implicit_solve(gamma,
                     reinterpret_cast<const CMFEM_Vector *>(&u),
                     reinterpret_cast<CMFEM_Vector *>(&k),
                     GetTime(),
                     implicit_context);
   }
};

CMFEM_ASSERT_TYPE(CMFEM_TimeDependentOperator, CTimeDependentOperatorAdapter);

mfem::TimeDependentOperator::Type AsMfemType(
   enum CMFEM_TimeDependentOperatorType type)
{
   switch (type)
   {
      case CMFEM_TIME_DEPENDENT_OPERATOR_EXPLICIT:
         return mfem::TimeDependentOperator::EXPLICIT;
      case CMFEM_TIME_DEPENDENT_OPERATOR_IMPLICIT:
         return mfem::TimeDependentOperator::IMPLICIT;
      case CMFEM_TIME_DEPENDENT_OPERATOR_HOMOGENEOUS:
         return mfem::TimeDependentOperator::HOMOGENEOUS;
   }
   return mfem::TimeDependentOperator::IMPLICIT;
}

mfem::TimeDependentOperator::ImplicitVariableType AsMfemType(
   enum CMFEM_ImplicitVariableType type)
{
   switch (type)
   {
      case CMFEM_IMPLICIT_VARIABLE_SLOPE:
         return mfem::TimeDependentOperator::ImplicitVariableType::SLOPE;
      case CMFEM_IMPLICIT_VARIABLE_STATE:
         return mfem::TimeDependentOperator::ImplicitVariableType::STATE;
   }
   return mfem::TimeDependentOperator::ImplicitVariableType::SLOPE;
}

} // namespace

extern "C" {

   CMFEM_TimeDependentOperator *CMFEM_TimeDependentOperator_New(
      int size,
      enum CMFEM_TimeDependentOperatorType type,
      CMFEM_TimeDependentOperatorMultCallback mult,
      void *mult_context,
      CMFEM_TimeDependentOperatorImplicitSolveCallback implicit_solve,
      void *implicit_context)
   {
      return reinterpret_cast<CMFEM_TimeDependentOperator *>(
                new CTimeDependentOperatorAdapter(size,
                                                  AsMfemType(type),
                                                  mult,
                                                  mult_context,
                                                  implicit_solve,
                                                  implicit_context));
   }

   void CMFEM_TimeDependentOperator_Delete(CMFEM_TimeDependentOperator *oper)
   {
      delete cmfem::As<CTimeDependentOperatorAdapter>(oper);
   }

   void CMFEM_TimeDependentOperator_SetImplicitVariableType(
      CMFEM_TimeDependentOperator *oper,
      enum CMFEM_ImplicitVariableType variable_type)
   {
      cmfem::As<mfem::TimeDependentOperator>(oper)->SetImplicitVariableType(
         AsMfemType(variable_type));
   }

} // extern "C"
