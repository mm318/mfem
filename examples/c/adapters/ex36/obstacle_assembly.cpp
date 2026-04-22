#include "common.hpp"

namespace
{

class LogObstacleCoefficient : public mfem::Coefficient
{
private:
   const mfem::GridFunction &state;
   mfem::Coefficient &obstacle;
   mfem::real_t min_value;

public:
   LogObstacleCoefficient(const mfem::GridFunction &state_,
                          mfem::Coefficient &obstacle_,
                          mfem::real_t min_value_)
      : state(state_), obstacle(obstacle_), min_value(min_value_) { }

   mfem::real_t Eval(mfem::ElementTransformation &T,
                     const mfem::IntegrationPoint &ip) override
   {
      const mfem::real_t value = state.GetValue(T, ip) - obstacle.Eval(T, ip);
      return std::max(min_value, std::log(value));
   }
};

class ExpObstacleCoefficient : public mfem::Coefficient
{
private:
   const mfem::GridFunction &state;
   mfem::Coefficient &obstacle;
   mfem::real_t min_value;
   mfem::real_t max_value;

public:
   ExpObstacleCoefficient(const mfem::GridFunction &state_,
                          mfem::Coefficient &obstacle_,
                          mfem::real_t min_value_,
                          mfem::real_t max_value_)
      : state(state_),
        obstacle(obstacle_),
        min_value(min_value_),
        max_value(max_value_) { }

   mfem::real_t Eval(mfem::ElementTransformation &T,
                     const mfem::IntegrationPoint &ip) override
   {
      const mfem::real_t value = std::exp(state.GetValue(T, ip))
                                 + obstacle.Eval(T, ip);
      return std::min(max_value, std::max(min_value, value));
   }
};

class NegExpCoefficient : public mfem::Coefficient
{
private:
   const mfem::GridFunction &state;
   mfem::real_t min_value;
   mfem::real_t max_value;

public:
   NegExpCoefficient(const mfem::GridFunction &state_,
                     mfem::real_t min_value_,
                     mfem::real_t max_value_)
      : state(state_), min_value(min_value_), max_value(max_value_) { }

   mfem::real_t Eval(mfem::ElementTransformation &T,
                     const mfem::IntegrationPoint &ip) override
   {
      const mfem::real_t value = std::exp(state.GetValue(T, ip));
      return -std::min(max_value, std::max(min_value, value));
   }
};

} // namespace

extern "C" {

   void CMFEM_ProjectObstacleLogGfFc(CMFEM_GridFunction *grid_function,
                                     const CMFEM_GridFunction *state,
                                     const CMFEM_FunctionCoefficient *obstacle,
                                     double min_val)
   {
      auto &obstacle_ref = const_cast<mfem::FunctionCoefficient &>(
                              *cmfem::As<const mfem::FunctionCoefficient>(
                                 obstacle));
      LogObstacleCoefficient coefficient(
         *cmfem::As<const mfem::GridFunction>(state),
         obstacle_ref,
         static_cast<mfem::real_t>(min_val));
      cmfem::As<mfem::GridFunction>(grid_function)->ProjectCoefficient(coefficient);
   }

   void CMFEM_ProjectObstacleExpGfFc(CMFEM_GridFunction *grid_function,
                                     const CMFEM_GridFunction *state,
                                     const CMFEM_FunctionCoefficient *obstacle,
                                     double min_val,
                                     double max_val)
   {
      auto &obstacle_ref = const_cast<mfem::FunctionCoefficient &>(
                              *cmfem::As<const mfem::FunctionCoefficient>(
                                 obstacle));
      ExpObstacleCoefficient coefficient(
         *cmfem::As<const mfem::GridFunction>(state),
         obstacle_ref,
         static_cast<mfem::real_t>(min_val),
         static_cast<mfem::real_t>(max_val));
      cmfem::As<mfem::GridFunction>(grid_function)->ProjectCoefficient(coefficient);
   }

   void CMFEM_AssembleObstacleLfGfFc(CMFEM_FiniteElementSpace *fespace,
                                     const CMFEM_GridFunction *state,
                                     const CMFEM_FunctionCoefficient *obstacle,
                                     CMFEM_Vector *out,
                                     double min_val,
                                     double max_val)
   {
      auto &obstacle_ref = const_cast<mfem::FunctionCoefficient &>(
                              *cmfem::As<const mfem::FunctionCoefficient>(
                                 obstacle));
      ExpObstacleCoefficient coefficient(
         *cmfem::As<const mfem::GridFunction>(state),
         obstacle_ref,
         static_cast<mfem::real_t>(min_val),
         static_cast<mfem::real_t>(max_val));
      mfem::LinearForm linear_form(cmfem::As<mfem::FiniteElementSpace>(fespace));
      linear_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(coefficient));
      linear_form.Assemble();
      cmfem::VectorRef(out) = linear_form;
   }

   CMFEM_SparseMatrix *CMFEM_AssembleObstacleHessianSmGf(
      CMFEM_FiniteElementSpace *fespace,
      const CMFEM_GridFunction *state,
      int order,
      double spectrum_shift,
      double min_val,
      double max_val)
   {
      NegExpCoefficient coefficient(*cmfem::As<const mfem::GridFunction>(state),
                                    static_cast<mfem::real_t>(min_val),
                                    static_cast<mfem::real_t>(max_val));
      mfem::ConstantCoefficient shift(static_cast<mfem::real_t>(spectrum_shift));
      mfem::BilinearForm bilinear_form(cmfem::As<mfem::FiniteElementSpace>(fespace));

      bilinear_form.AddDomainIntegrator(new mfem::MassIntegrator(coefficient));
      if (order == 1)
      {
         bilinear_form.AddDomainIntegrator(new mfem::MassIntegrator(shift));
      }
      else
      {
         bilinear_form.AddDomainIntegrator(new mfem::DiffusionIntegrator(shift));
      }
      bilinear_form.Assemble();
      bilinear_form.Finalize();

      return reinterpret_cast<CMFEM_SparseMatrix *>(
                new mfem::SparseMatrix(bilinear_form.SpMat()));
   }

} // extern "C"
