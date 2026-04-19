#include "common.hpp"

#include <fstream>
#include <sstream>

namespace
{

void ComputeCurrentDensityOnSubMesh(mfem::Mesh &mesh_cond,
                                    int order,
                                    bool pa,
                                    const mfem::Array<int> &phi0_attr,
                                    const mfem::Array<int> &phi1_attr,
                                    const mfem::Array<int> &jn_zero_attr,
                                    mfem::GridFunction &j_cond)
{
   mfem::FiniteElementSpace &fes_cond_rt = *j_cond.FESpace();
   const int dim = mesh_cond.Dimension();

   mfem::H1_FECollection fec_h1(order, dim);
   mfem::FiniteElementSpace fes_cond_h1(&mesh_cond, &fec_h1);

   mfem::ConstantCoefficient sigma_coef(1.0);
   mfem::Array<int> ess_bdr_phi(mesh_cond.bdr_attributes.Max());
   mfem::Array<int> ess_bdr_j(mesh_cond.bdr_attributes.Max());
   mfem::Array<int> ess_bdr_tdof_phi;
   ess_bdr_phi = 0;
   ess_bdr_j = 0;
   for (int i = 0; i < phi0_attr.Size(); i++)
   {
      ess_bdr_phi[phi0_attr[i] - 1] = 1;
   }
   for (int i = 0; i < phi1_attr.Size(); i++)
   {
      ess_bdr_phi[phi1_attr[i] - 1] = 1;
   }
   for (int i = 0; i < jn_zero_attr.Size(); i++)
   {
      ess_bdr_j[jn_zero_attr[i] - 1] = 1;
   }
   fes_cond_h1.GetEssentialTrueDofs(ess_bdr_phi, ess_bdr_tdof_phi);

   mfem::BilinearForm a_h1(&fes_cond_h1);
   a_h1.AddDomainIntegrator(new mfem::DiffusionIntegrator(sigma_coef));
   a_h1.Assemble();

   mfem::LinearForm b_h1(&fes_cond_h1);
   b_h1 = 0.0;

   mfem::ConstantCoefficient one(1.0);
   mfem::ConstantCoefficient zero(0.0);
   mfem::GridFunction phi_h1(&fes_cond_h1);
   phi_h1 = 0.0;

   mfem::Array<int> bdr0(mesh_cond.bdr_attributes.Max());
   bdr0 = 0;
   for (int i = 0; i < phi0_attr.Size(); i++)
   {
      bdr0[phi0_attr[i] - 1] = 1;
   }
   phi_h1.ProjectBdrCoefficient(zero, bdr0);

   mfem::Array<int> bdr1(mesh_cond.bdr_attributes.Max());
   bdr1 = 0;
   for (int i = 0; i < phi1_attr.Size(); i++)
   {
      bdr1[phi1_attr[i] - 1] = 1;
   }
   phi_h1.ProjectBdrCoefficient(one, bdr1);

   {
      mfem::OperatorPtr A;
      mfem::Vector B, X;
      a_h1.FormLinearSystem(ess_bdr_tdof_phi, phi_h1, b_h1, A, X, B);

      if (!pa)
      {
         mfem::out << "\nSolving for electric potential using PCG "
                   << "with a Gauss-Seidel preconditioner" << std::endl;
         mfem::GSSmoother M((mfem::SparseMatrix &)(*A));
         mfem::PCG(*A, M, B, X, 1, 200, 1e-12, 0.0);
      }
      else
      {
         mfem::out << "\nSolving for electric potential using CG" << std::endl;
         if (mfem::UsesTensorBasis(fes_cond_h1))
         {
            mfem::OperatorJacobiSmoother M(a_h1, ess_bdr_tdof_phi);
            mfem::PCG(*A, M, B, X, 1, 400, 1e-12, 0.0);
         }
         else
         {
            mfem::CG(*A, B, X, 1, 400, 1e-12, 0.0);
         }
      }
      a_h1.RecoverFEMSolution(X, b_h1, phi_h1);
   }

   mfem::BilinearForm m_rt(&fes_cond_rt);
   m_rt.AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
   m_rt.Assemble();

   mfem::MixedBilinearForm d_h1(&fes_cond_h1, &fes_cond_rt);
   d_h1.AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(sigma_coef));
   d_h1.Assemble();

   mfem::LinearForm b_rt(&fes_cond_rt);
   d_h1.Mult(phi_h1, b_rt);
   b_rt *= -1.0;

   mfem::out << "\nSolving for current density in H(Div) "
             << "using diagonally scaled CG" << std::endl;
   mfem::out << "Size of linear system: "
             << fes_cond_rt.GetTrueVSize() << std::endl;

   mfem::Array<int> ess_bdr_tdof_rt;
   mfem::OperatorPtr M;
   mfem::Vector B, X;

   fes_cond_rt.GetEssentialTrueDofs(ess_bdr_j, ess_bdr_tdof_rt);

   j_cond = 0.0;
   m_rt.FormLinearSystem(ess_bdr_tdof_rt, j_cond, b_rt, M, X, B);

   mfem::CGSolver cg;
   cg.SetRelTol(1e-12);
   cg.SetMaxIter(2000);
   cg.SetPrintLevel(1);
   cg.SetOperator(*M);
   cg.Mult(B, X);
   m_rt.RecoverFEMSolution(X, b_rt, j_cond);
}

} // namespace

extern "C" {

   void CMFEM_Ex34_PrepareCurrentDensityRt(CMFEM_Mesh *mesh,
                                           int order,
                                           int pa,
                                           int visualization,
                                           const CMFEM_ArrayInt *submesh_elems,
                                           const CMFEM_ArrayInt *phi0_attr,
                                           const CMFEM_ArrayInt *phi1_attr,
                                           const CMFEM_ArrayInt *jn_zero_attr,
                                           CMFEM_GridFunction *j_full)
   {
      mfem::Mesh &mesh_ref = *cmfem::As<mfem::Mesh>(mesh);
      mfem::Array<int> cond_attr(1);
      const int submesh_attr = mesh_ref.attributes.Max() + 1;
      cond_attr[0] = submesh_attr;

      for (int i = 0; i < cmfem::ArrayIntRef(submesh_elems).Size(); i++)
      {
         mesh_ref.SetAttribute(cmfem::ArrayIntRef(submesh_elems)[i], submesh_attr);
      }
      mesh_ref.SetAttributes();

      mfem::SubMesh mesh_cond(mfem::SubMesh::CreateFromDomain(mesh_ref, cond_attr));

      mfem::RT_FECollection fec_cond_rt(order - 1, mesh_cond.Dimension());
      mfem::FiniteElementSpace fes_cond_rt(&mesh_cond, &fec_cond_rt);
      mfem::GridFunction j_cond(&fes_cond_rt);

      ComputeCurrentDensityOnSubMesh(mesh_cond,
                                     order,
                                     pa != 0,
                                     cmfem::ArrayIntRef(phi0_attr),
                                     cmfem::ArrayIntRef(phi1_attr),
                                     cmfem::ArrayIntRef(jn_zero_attr),
                                     j_cond);

      {
         std::ofstream mesh_ofs("cond.mesh");
         mesh_ofs.precision(8);
         mesh_cond.Print(mesh_ofs);

         std::ofstream cond_ofs("cond_j.gf");
         cond_ofs.precision(8);
         j_cond.Save(cond_ofs);
      }

      if (visualization)
      {
         char vishost[] = "localhost";
         int visport = 19916;
         mfem::socketstream port_sock(vishost, visport);
         port_sock.precision(8);
         port_sock << "solution\n" << mesh_cond << j_cond << std::flush;
      }

      *cmfem::As<mfem::GridFunction>(j_full) = 0.0;
      mesh_cond.Transfer(j_cond, *cmfem::As<mfem::GridFunction>(j_full));
   }

   void CMFEM_Ex34_AssembleCurrentDensityLf(CMFEM_FiniteElementSpace *fespace_nd,
                                            const CMFEM_GridFunction *j_full,
                                            CMFEM_LinearForm *linear_form)
   {
      mfem::VectorGridFunctionCoefficient coef(
         const_cast<mfem::GridFunction *>(cmfem::As<const mfem::GridFunction>(
                                             j_full)));
      auto *lf = cmfem::As<mfem::LinearForm>(linear_form);
      (void)fespace_nd;
      lf->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(coef));
      lf->Assemble();
   }

} // extern "C"
