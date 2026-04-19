#include "common.hpp"

#include "config/runtime_paths.hpp"
#include "examples/cpp/ex18.hpp"

#include <fstream>
#include <memory>
#include <sstream>

namespace
{

class Ex18DghclContext
{
private:
   mfem::EulerFlux flux_;
   mfem::RusanovFlux numerical_flux_;
   mfem::DGHyperbolicConservationLaws oper_;

public:
   Ex18DghclContext(mfem::FiniteElementSpace &vfes,
                    int int_order_offset,
                    mfem::real_t specific_heat_ratio,
                    bool preassemble_weak_div)
      : flux_(vfes.GetMesh()->Dimension(), specific_heat_ratio),
        numerical_flux_(flux_),
        oper_(
           vfes,
           std::unique_ptr<mfem::HyperbolicFormIntegrator>(
              new mfem::HyperbolicFormIntegrator(numerical_flux_,
                                                 int_order_offset)),
           preassemble_weak_div)
   { }

   mfem::DGHyperbolicConservationLaws &Operator()
   {
      return oper_;
   }
};

void SaveStateImpl(const mfem::Mesh &mesh,
                   const mfem::FiniteElementSpace &scalar_fes,
                   int num_equations,
                   const mfem::GridFunction &sol,
                   const char *mesh_path,
                   const char *tag)
{
   std::ofstream mesh_ofs(mesh_path);
   mesh_ofs.precision(8);
   mesh.Print(mesh_ofs);

   const int scalar_ndofs = scalar_fes.GetNDofs();
   for (int k = 0; k < num_equations; k++)
   {
      mfem::GridFunction uk(const_cast<mfem::FiniteElementSpace *>(&scalar_fes),
                            const_cast<mfem::real_t *>(sol.GetData()) +
                            k * scalar_ndofs);
      std::ostringstream sol_name;
      sol_name << "euler-" << k << "-" << tag << ".gf";
      std::ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      uk.Save(sol_ofs);
   }
}

} // namespace

extern "C" {

   CMFEM_Mesh *CMFEM_Ex18_NewDefaultMesh(int problem)
   {
      return reinterpret_cast<CMFEM_Mesh *>(new mfem::Mesh(mfem::EulerMesh(problem)));
   }

   void CMFEM_Ex18_ProjectInitialConditionGf(int problem,
                                             double specific_heat_ratio,
                                             double gas_constant,
                                             CMFEM_GridFunction *sol)
   {
      auto u0 = mfem::EulerInitialCondition(problem,
                                            static_cast<mfem::real_t>(
                                               specific_heat_ratio),
                                            static_cast<mfem::real_t>(
                                               gas_constant));
      cmfem::As<mfem::GridFunction>(sol)->ProjectCoefficient(u0);
   }

   void *CMFEM_Ex18_NewDghcl(CMFEM_FiniteElementSpace *vfes,
                             int int_order_offset,
                             double specific_heat_ratio,
                             int preassemble_weak_div)
   {
      return new Ex18DghclContext(
                *cmfem::As<mfem::FiniteElementSpace>(vfes),
                int_order_offset,
                static_cast<mfem::real_t>(specific_heat_ratio),
                preassemble_weak_div != 0);
   }

   void CMFEM_Ex18_DeleteDghcl(void *context)
   {
      delete static_cast<Ex18DghclContext *>(context);
   }

   void CMFEM_Ex18_SetTimeDghcl(void *context, double time)
   {
      static_cast<Ex18DghclContext *>(context)->Operator().SetTime(
         static_cast<mfem::real_t>(time));
   }

   void CMFEM_Ex18_MultDghcl(void *context,
                             const CMFEM_Vector *x,
                             CMFEM_Vector *y)
   {
      static_cast<Ex18DghclContext *>(context)->Operator().Mult(
         cmfem::VectorRef(x),
         cmfem::VectorRef(y));
   }

   double CMFEM_Ex18_GetMaxCharSpeedDghcl(void *context)
   {
      return static_cast<double>(
                static_cast<Ex18DghclContext *>(context)->Operator()
                .GetMaxCharSpeed());
   }

   double CMFEM_Ex18_ComputeStableDtMeshDghclGf(const CMFEM_Mesh *mesh,
                                                int order,
                                                double cfl,
                                                void *context,
                                                const CMFEM_GridFunction *sol)
   {
      mfem::real_t hmin = mfem::infinity();
      mfem::Vector x;
      mfem::Vector z(cmfem::As<const mfem::GridFunction>(sol)->Size());
      mfem::Mesh &mesh_ref = *const_cast<mfem::Mesh *>(
                                cmfem::As<const mfem::Mesh>(mesh));

      for (int i = 0; i < mesh_ref.GetNE(); i++)
      {
         hmin = std::min(mesh_ref.GetElementSize(i, 1), hmin);
      }

      x = *cmfem::As<const mfem::GridFunction>(sol);
      static_cast<Ex18DghclContext *>(context)->Operator().Mult(x, z);

      return static_cast<double>(
                static_cast<mfem::real_t>(cfl) *
                hmin /
                static_cast<Ex18DghclContext *>(context)->Operator()
                .GetMaxCharSpeed() /
                static_cast<mfem::real_t>(2 * order + 1));
   }

   void CMFEM_Ex18_SaveStateMeshSfesGf(const CMFEM_Mesh *mesh,
                                       const CMFEM_FiniteElementSpace *scalar_fes,
                                       int num_equations,
                                       const CMFEM_GridFunction *sol,
                                       const char *mesh_path,
                                       const char *tag)
   {
      SaveStateImpl(*cmfem::As<const mfem::Mesh>(mesh),
                    *cmfem::As<const mfem::FiniteElementSpace>(scalar_fes),
                    num_equations,
                    *cmfem::As<const mfem::GridFunction>(sol),
                    mesh_path,
                    tag);
   }

   double CMFEM_Ex18_ComputeL2ErrorGf(int problem,
                                      double specific_heat_ratio,
                                      double gas_constant,
                                      const CMFEM_GridFunction *sol)
   {
      auto u0 = mfem::EulerInitialCondition(problem,
                                            static_cast<mfem::real_t>(
                                               specific_heat_ratio),
                                            static_cast<mfem::real_t>(
                                               gas_constant));
      return static_cast<double>(
                cmfem::As<const mfem::GridFunction>(sol)->ComputeLpError(2, u0));
   }

   int CMFEM_Ex18_SendMomentumToGLVisMeshSfesDfesGf(
      const CMFEM_Mesh *mesh,
      const CMFEM_FiniteElementSpace *scalar_fes,
      const CMFEM_FiniteElementSpace *dfes,
      const CMFEM_GridFunction *sol,
      double time,
      int pause)
   {
      const int scalar_ndofs =
         cmfem::As<const mfem::FiniteElementSpace>(scalar_fes)->GetNDofs();
      mfem::GridFunction mom(
         const_cast<mfem::FiniteElementSpace *>(cmfem::As<
                                                const mfem::FiniteElementSpace>(dfes)),
         const_cast<mfem::real_t *>(cmfem::As<const mfem::GridFunction>(sol)
                                    ->GetData()) +
         scalar_ndofs);

      char vishost[] = "localhost";
      int visport = 19916;
      mfem::socketstream sout(vishost, visport);
      if (!sout)
      {
         return 0;
      }

      sout.precision(8);
      sout << "solution\n"
           << *cmfem::As<const mfem::Mesh>(mesh)
           << mom;
      sout << "window_title 'momentum, t = " << time << "'\n";
      sout << "view 0 0\n";
      sout << "keys jlm\n";
      if (pause)
      {
         sout << "pause\n";
      }
      sout << std::flush;
      return 1;
   }

} // extern "C"
