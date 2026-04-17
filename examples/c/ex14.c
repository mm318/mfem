//                              MFEM C Example 14
//
// Description: This example code demonstrates a discontinuous Galerkin (DG)
//              finite element discretization of the Poisson problem
//              -Delta u = 1 with homogeneous Dirichlet boundary conditions.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int ref_levels = -1;
   int order = 1;
   double sigma = -1.0;
   double kappa = -1.0;
   double eta = 0.0;
   int pa = 0;
   int visualization = 1;
   const char *device_config = "cpu";
   int dim;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-r", "--refine",
                                      &ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

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

      parsed = cmfem_parse_double_option(argc, argv, &i, "-e", "--eta", &eta);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-pa", "--partial-assembly",
                                       "-no-pa", "--no-partial-assembly",
                                       &pa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device",
                                         &device_config);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   if (kappa < 0.0)
   {
      kappa = (double)(order + 1) * (double)(order + 1);
   }

   // 2. Enable hardware devices such as GPUs and programming models such as
   //    OpenMP based on command line options.
   CMFEM_Device *device = CMFEM_Device_New(device_config);
   CMFEM_Device_Print(device);

   // 3. Read the mesh from the given mesh file.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   dim = CMFEM_Mesh_Dimension(mesh);

   // 4. Refine the mesh to increase the resolution and project NURBS meshes to
   //    piecewise-polynomial curved meshes.
   if (ref_levels < 0)
   {
      ref_levels = cmfem_uniform_refinement_levels(50000.0, CMFEM_Mesh_GetNE(mesh),
                                                   dim);
   }
   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }
   if (CMFEM_Mesh_HasNURBSext(mesh))
   {
      CMFEM_Mesh_SetCurvature(mesh, order > 0 ? order : 1);
   }

   // 5. Define a discontinuous finite element space on the mesh.
   CMFEM_DG_FECollection *fec =
      CMFEM_DG_FECollection_NewOrderDimBasis(order, dim,
                                             pa ? CMFEM_BASIS_GAUSS_LOBATTO :
                                             CMFEM_BASIS_GAUSS_LEGENDRE);
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshDG(mesh,
                                                                          fec);
   printf("Number of unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   // 6. Set up the linear form b(.) for the right-hand side and weakly imposed
   //    homogeneous Dirichlet boundary condition.
   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_ConstantCoefficient *zero = CMFEM_ConstantCoefficient_New(0.0);
   CMFEM_LinearForm_AddDomainIntegrator_DomainLFIntegrator_ConstantCoefficient(b,
                                                                               one);
   CMFEM_LinearForm_AddBdrFaceIntegrator_DGDirichletLFIntegrator(
      b, zero, one, sigma, kappa);
   CMFEM_LinearForm_Assemble(b);

   // 7. Define the solution vector x corresponding to fespace and initialize it
   //    with zero.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_Assign(x, 0.0);

   // 8. Set up the bilinear form for the DG discretization.
   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   CMFEM_BilinearForm_AddDomainIntegrator_DiffusionCoefficient(a, one);
   CMFEM_BilinearForm_AddInteriorFaceIntegrator_DGDiffusionIntegrator(
      a, one, sigma, kappa);
   CMFEM_BilinearForm_AddBdrFaceIntegrator_DGDiffusionIntegrator(
      a, one, sigma, kappa);
   if (eta > 0.0)
   {
      if (pa)
      {
         fprintf(stderr, "BR2 is not compatible with partial assembly.\n");
         CMFEM_BilinearForm_Delete(a);
         CMFEM_GridFunction_Delete(x);
         CMFEM_ConstantCoefficient_Delete(zero);
         CMFEM_ConstantCoefficient_Delete(one);
         CMFEM_LinearForm_Delete(b);
         CMFEM_FiniteElementSpace_Delete(fespace);
         CMFEM_DG_FECollection_Delete(fec);
         CMFEM_Mesh_Delete(mesh);
         CMFEM_Device_Delete(device);
         return 1;
      }
      CMFEM_BilinearForm_AddInteriorFaceIntegrator_DGDiffusionBR2Integrator(
         a, fespace, eta);
      CMFEM_BilinearForm_AddBdrFaceIntegrator_DGDiffusionBR2Integrator(
         a, fespace, eta);
   }
   if (pa)
   {
      CMFEM_BilinearForm_SetAssemblyLevelPartial(a);
   }
   CMFEM_BilinearForm_Assemble(a);
   CMFEM_BilinearForm_Finalize(a);

   // 9. Assemble the linear system and solve it with CG or GMRES depending on
   //    the symmetry of the DG operator.
   _Alignas(max_align_t) CMFEM_ArrayInt empty_tdof_list =
      CMFEM_ArrayInt_Construct();
   _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
   CMFEM_BilinearForm_FormLinearSystemOperator(a, &empty_tdof_list, x, b, &A, &X,
                                               &B);
   if (pa && sigma != -1.0)
   {
      fprintf(stderr,
              "Partial assembly currently requires sigma = -1 in this C port.\n");
      CMFEM_OperatorPtr_Destroy(&A);
      CMFEM_Vector_Destroy(&B);
      CMFEM_Vector_Destroy(&X);
      CMFEM_ArrayInt_Destroy(&empty_tdof_list);
      CMFEM_BilinearForm_Delete(a);
      CMFEM_GridFunction_Delete(x);
      CMFEM_ConstantCoefficient_Delete(zero);
      CMFEM_ConstantCoefficient_Delete(one);
      CMFEM_LinearForm_Delete(b);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_DG_FECollection_Delete(fec);
      CMFEM_Mesh_Delete(mesh);
      CMFEM_Device_Delete(device);
      return 1;
   }

   if (pa)
   {
      CMFEM_CG_Operator(&A, &B, &X, 1, 500, 1e-12, 0.0);
   }
   else
   {
      CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewOperator(&A);
      if (sigma == -1.0)
      {
         CMFEM_PCG_OperatorGSSmoother(&A, M, &B, &X, 1, 500, 1e-12, 0.0);
      }
      else
      {
         CMFEM_GMRES_OperatorGSSmoother(&A, M, &B, &X, 1, 500, 10, 1e-12, 0.0);
      }
      CMFEM_GSSmoother_Delete(M);
   }

   // 10. Recover the solution as a finite element grid function.
   CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);

   // 11. Save the refined mesh and the solution.
   CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
   CMFEM_GridFunction_Save(x, "sol.gf", 8);

   // 12. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
   }

   CMFEM_OperatorPtr_Destroy(&A);
   CMFEM_Vector_Destroy(&B);
   CMFEM_Vector_Destroy(&X);
   CMFEM_ArrayInt_Destroy(&empty_tdof_list);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_GridFunction_Delete(x);
   CMFEM_ConstantCoefficient_Delete(zero);
   CMFEM_ConstantCoefficient_Delete(one);
   CMFEM_LinearForm_Delete(b);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_DG_FECollection_Delete(fec);
   CMFEM_Mesh_Delete(mesh);
   CMFEM_Device_Delete(device);
   return 0;
}
