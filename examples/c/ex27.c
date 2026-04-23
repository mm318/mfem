//                              MFEM C Example 27
//
// Description: This example solves a Laplace problem with Dirichlet, Neumann,
//              Robin, and periodic boundary conditions using either a
//              continuous H1 basis or a discontinuous Galerkin basis.

#include "common.h"
#include "cmfem.h"
#include "adapters/ex27/boundary_postprocess.h"
#include "adapters/ex27/periodic_two_hole_mesh.h"

#include <stddef.h>
#include <stdio.h>

static void report_boundary_condition(const CMFEM_GridFunction *solution,
                                      const CMFEM_ArrayInt *marker,
                                      const char *label,
                                      double alpha,
                                      double beta,
                                      double gamma)
{
   double error = 0.0;
   double avg = CMFEM_IntegrateBoundaryConditionGfAi(solution,
                                                     marker,
                                                     alpha,
                                                     beta,
                                                     gamma,
                                                     &error);
   const int homogeneous = gamma == 0.0;
   if (!homogeneous)
   {
      error /= fabs(gamma);
   }

   printf("Average of %s:\t%g, \t%s error %g\n",
          label,
          avg,
          homogeneous ? "absolute" : "relative",
          error);
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int ser_ref_levels = 2;
   int order = 1;
   double sigma = -1.0;
   double kappa = -1.0;
   int h1 = 1;
   int visualization = 1;
   double mat_val = 1.0;
   double dbc_val = 0.0;
   double nbc_val = 1.0;
   double rbc_a_val = 1.0;
   double rbc_b_val = 1.0;
   double radius = 0.2;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_bool_option(&i, argv,
                                           "-h1", "--continuous",
                                           "-dg", "--discontinuous",
                                           &h1);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-s", "--sigma",
                                         &sigma);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-k", "--kappa",
                                         &kappa);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-rs",
                                      "--refine-serial", &ser_ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-mat",
                                         "--material-value", &mat_val);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-dbc",
                                         "--dirichlet-value", &dbc_val);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-nbc",
                                         "--neumann-value", &nbc_val);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-rbc-a",
                                         "--robin-a-value", &rbc_a_val);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-rbc-b",
                                         "--robin-b-value", &rbc_b_val);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-a", "--radius",
                                         &radius);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   if (kappa < 0.0 && !h1)
   {
      kappa = (double)((order + 1) * (order + 1));
   }
   if (radius < 0.01)
   {
      printf("Hole radius too small, resetting to 0.01.\n");
      radius = 0.01;
   }
   if (radius > 0.49)
   {
      printf("Hole radius too large, resetting to 0.49.\n");
      radius = 0.49;
   }

   // 2. Construct the periodic two-hole mesh and the finite element space.
   CMFEM_Mesh *mesh = CMFEM_NewPeriodicTwoHoleMesh2d(ser_ref_levels, radius);
   CMFEM_H1FeCollection *h1_fec = NULL;
   CMFEM_DgFeCollection *dg_fec = NULL;
   CMFEM_FiniteElementSpace *fespace = NULL;
   CMFEM_GridFunction *u = NULL;
   CMFEM_BilinearForm *a = NULL;
   CMFEM_LinearForm *b = NULL;
   CMFEM_ConstantCoefficient *mat_coef = NULL;
   CMFEM_ConstantCoefficient *dbc_coef = NULL;
   CMFEM_ConstantCoefficient *m_nbc_coef = NULL;
   CMFEM_ConstantCoefficient *m_rbc_a_coef = NULL;
   CMFEM_ConstantCoefficient *m_rbc_b_coef = NULL;
   CMFEM_GSSmoother *smoother = NULL;
   _Alignas(max_align_t) CMFEM_ArrayInt nbc_bdr =
      CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
   _Alignas(max_align_t) CMFEM_ArrayInt rbc_bdr =
      CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
   _Alignas(max_align_t) CMFEM_ArrayInt dbc_bdr =
      CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
   _Alignas(max_align_t) CMFEM_ArrayInt nbc0_bdr =
      CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
   _Alignas(max_align_t) CMFEM_SparseMatrix A = CMFEM_SparseMatrix_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();

   if (h1)
   {
      h1_fec = CMFEM_H1FeCollection_NewOrderDim(order, CMFEM_Mesh_Dimension(mesh));
      fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, h1_fec);
   }
   else
   {
      dg_fec = CMFEM_DgFeCollection_NewOrderDimBasis(order,
                                                     CMFEM_Mesh_Dimension(mesh),
                                                     CMFEM_BASIS_GAUSS_LEGENDRE);
      fespace = CMFEM_FiniteElementSpace_NewMeshDg(mesh, dg_fec);
   }
   printf("Number of finite element unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   // 3. Mark the boundary attributes associated with each boundary condition.
   CMFEM_ArrayInt_Assign(&nbc_bdr, 0);
   CMFEM_ArrayInt_Assign(&rbc_bdr, 0);
   CMFEM_ArrayInt_Assign(&dbc_bdr, 0);
   CMFEM_ArrayInt_Assign(&nbc0_bdr, 0);
   CMFEM_ArrayInt_Set(&nbc_bdr, 0, 1);
   CMFEM_ArrayInt_Set(&rbc_bdr, 1, 1);
   CMFEM_ArrayInt_Set(&dbc_bdr, 2, 1);
   CMFEM_ArrayInt_Set(&nbc0_bdr, 3, 1);
   if (h1 && CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
   {
      CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace,
                                                    &dbc_bdr,
                                                    &ess_tdof_list);
   }

   // 4. Set up the coefficients, solution grid function, and discrete forms.
   mat_coef = CMFEM_ConstantCoefficient_New(mat_val);
   dbc_coef = CMFEM_ConstantCoefficient_New(dbc_val);
   m_nbc_coef = CMFEM_ConstantCoefficient_New(mat_val * nbc_val);
   m_rbc_a_coef = CMFEM_ConstantCoefficient_New(mat_val * rbc_a_val);
   m_rbc_b_coef = CMFEM_ConstantCoefficient_New(mat_val * rbc_b_val);

   u = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_Assign(u, 0.0);

   a = CMFEM_BilinearForm_New(fespace);
   CMFEM_BilinearForm_AddDomainIntegratorDiCc(a, mat_coef);
   if (h1)
   {
      CMFEM_BilinearForm_AddBoundaryIntegratorMiCcAi(a, m_rbc_a_coef, &rbc_bdr);
   }
   else
   {
      CMFEM_BilinearForm_AddInteriorFaceIntegratorDgd(a, mat_coef, sigma, kappa);
      CMFEM_BilinearForm_AddBdrFaceIntegratorDgdAi(a, mat_coef, sigma, kappa,
                                                   &dbc_bdr);
      CMFEM_BilinearForm_AddBdrFaceIntegratorBmiCcAi(a, m_rbc_a_coef,
                                                     &rbc_bdr);
   }
   CMFEM_BilinearForm_Assemble(a);
   CMFEM_BilinearForm_Finalize(a);

   b = CMFEM_LinearForm_New(fespace);
   if (h1)
   {
      CMFEM_GridFunction_ProjectBdrCoefficientCc(u, dbc_coef, &dbc_bdr);
      CMFEM_LinearForm_AddBoundaryIntegratorBliCcAi(b, m_nbc_coef, &nbc_bdr);
      CMFEM_LinearForm_AddBoundaryIntegratorBliCcAi(b, m_rbc_b_coef, &rbc_bdr);
   }
   else
   {
      CMFEM_LinearForm_AddBdrFaceIntegratorDglAi(b, dbc_coef, mat_coef,
                                                 sigma, kappa, &dbc_bdr);
      CMFEM_LinearForm_AddBdrFaceIntegratorBliCcAi(b, m_nbc_coef, &nbc_bdr);
      CMFEM_LinearForm_AddBdrFaceIntegratorBliCcAi(b, m_rbc_b_coef, &rbc_bdr);
   }
   CMFEM_LinearForm_Assemble(b);

   // 5. Form the linear system, solve it, and recover the finite element
   //    solution.
   CMFEM_BilinearForm_FormLinearSystemSm(a, &ess_tdof_list, u, b, &A, &X, &B);
   smoother = CMFEM_GSSmoother_NewSm(&A);
   if (sigma == -1.0)
   {
      CMFEM_PCGSmGs(&A, smoother, &B, &X, 1, 500, 1e-12, 0.0);
   }
   else
   {
      CMFEM_GMRESSmGs(&A, smoother, &B, &X, 1, 500, 10, 1e-12, 0.0);
   }
   CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, u);

   // 6. Verify the imposed boundary conditions by integrating the recovered
   //    solution and flux traces over the marked boundaries.
   printf("\nVerifying boundary conditions\n");
   printf("=============================\n");
   report_boundary_condition(u, &dbc_bdr, "solution on Gamma_dbc",
                             0.0, 1.0, dbc_val);
   report_boundary_condition(u, &nbc_bdr, "n.Grad(u) on Gamma_nbc",
                             1.0, 0.0, nbc_val);
   report_boundary_condition(u, &nbc0_bdr, "n.Grad(u) on Gamma_nbc0",
                             1.0, 0.0, 0.0);
   report_boundary_condition(u, &rbc_bdr, "n.Grad(u)+a*u on Gamma_rbc",
                             1.0, rbc_a_val, rbc_b_val);

   // 7. Save the refined mesh and the solution for offline visualization.
   CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
   CMFEM_GridFunction_Save(u, "sol.gf", 8);

   // 8. Optionally send the solution to GLVis.
   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, u, "localhost", 19916);
   }

   // 9. Free the allocated resources.
   CMFEM_GSSmoother_Delete(smoother);
   CMFEM_Vector_Destroy(&X);
   CMFEM_Vector_Destroy(&B);
   CMFEM_SparseMatrix_Destroy(&A);
   CMFEM_ArrayInt_Destroy(&ess_tdof_list);
   CMFEM_ArrayInt_Destroy(&nbc0_bdr);
   CMFEM_ArrayInt_Destroy(&dbc_bdr);
   CMFEM_ArrayInt_Destroy(&rbc_bdr);
   CMFEM_ArrayInt_Destroy(&nbc_bdr);
   CMFEM_ConstantCoefficient_Delete(m_rbc_b_coef);
   CMFEM_ConstantCoefficient_Delete(m_rbc_a_coef);
   CMFEM_ConstantCoefficient_Delete(m_nbc_coef);
   CMFEM_ConstantCoefficient_Delete(dbc_coef);
   CMFEM_ConstantCoefficient_Delete(mat_coef);
   CMFEM_LinearForm_Delete(b);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_GridFunction_Delete(u);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_DgFeCollection_Delete(dg_fec);
   CMFEM_H1FeCollection_Delete(h1_fec);
   CMFEM_Mesh_Delete(mesh);

   return 0;
}
