#include "common.hpp"

namespace
{

class CSecondOrderTimeDependentOperator final
   : public mfem::SecondOrderTimeDependentOperator
{
private:
   CMFEM_SecondOrderTimeDependentOperatorMultCallback mult;
   void *mult_context;
   CMFEM_SecondOrderTimeDependentOperatorImplicitSolveCallback implicit_solve;
   void *implicit_context;

public:
   CSecondOrderTimeDependentOperator(
      int size,
      CMFEM_SecondOrderTimeDependentOperatorMultCallback mult_,
      void *mult_context_,
      CMFEM_SecondOrderTimeDependentOperatorImplicitSolveCallback
      implicit_solve_,
      void *implicit_context_)
      : mfem::SecondOrderTimeDependentOperator(size, 0.0),
        mult(mult_),
        mult_context(mult_context_),
        implicit_solve(implicit_solve_),
        implicit_context(implicit_context_) { }

   void Mult(const mfem::Vector &u,
             const mfem::Vector &du_dt,
             mfem::Vector &d2udt2) const override
   {
      if (!mult)
      {
         MFEM_ABORT("CMFEM second-order operator missing Mult callback");
      }
      mult(reinterpret_cast<const CMFEM_Vector *>(&u),
           reinterpret_cast<const CMFEM_Vector *>(&du_dt),
           reinterpret_cast<CMFEM_Vector *>(&d2udt2),
           GetTime(),
           mult_context);
   }

   void ImplicitSolve(const mfem::real_t fac0,
                      const mfem::real_t fac1,
                      const mfem::Vector &u,
                      const mfem::Vector &du_dt,
                      mfem::Vector &d2udt2) override
   {
      if (!implicit_solve)
      {
         MFEM_ABORT(
            "CMFEM second-order operator missing ImplicitSolve callback");
      }
      implicit_solve(fac0,
                     fac1,
                     reinterpret_cast<const CMFEM_Vector *>(&u),
                     reinterpret_cast<const CMFEM_Vector *>(&du_dt),
                     reinterpret_cast<CMFEM_Vector *>(&d2udt2),
                     GetTime(),
                     implicit_context);
   }
};

CMFEM_ASSERT_TYPE(CMFEM_SecondOrderTimeDependentOperator,
                  CSecondOrderTimeDependentOperator);

} // namespace

extern "C" {

   CMFEM_SecondOrderTimeDependentOperator *
   CMFEM_SecondOrderTimeDependentOperator_New(
      int size,
      CMFEM_SecondOrderTimeDependentOperatorMultCallback mult,
      void *mult_context,
      CMFEM_SecondOrderTimeDependentOperatorImplicitSolveCallback
      implicit_solve,
      void *implicit_context)
   {
      return reinterpret_cast<CMFEM_SecondOrderTimeDependentOperator *>(
                new CSecondOrderTimeDependentOperator(size,
                                                      mult,
                                                      mult_context,
                                                      implicit_solve,
                                                      implicit_context));
   }

   void CMFEM_SecondOrderTimeDependentOperator_Delete(
      CMFEM_SecondOrderTimeDependentOperator *oper)
   {
      delete cmfem::As<CSecondOrderTimeDependentOperator>(oper);
   }

} // extern "C"
