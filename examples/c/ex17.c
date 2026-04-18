//                              MFEM C Example 17
//
// Description: This example solves a linear elasticity cantilever beam using
//              a discontinuous Galerkin vector finite element discretization.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

static void init_displacement(const CMFEM_Vector *x, CMFEM_Vector *u,
                              void *context)
{
   int i;
   (void)context;
   for (i = 0; i < CMFEM_Vector_Size(u); i++)
   {
      CMFEM_Vector_Set(u, i, 0.0);
   }
   CMFEM_Vector_Set(u, CMFEM_Vector_Size(u) - 1, -0.2 * CMFEM_Vector_Get(x, 0));
}

int main(int argc, char *argv[])
{
   // 1. Define and parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("beam-tri.mesh");
   int ref_levels = -1;
   int order = 1;
   double alpha = -1.0;
   double kappa = -1.0;
   int visualization = 1;
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

      parsed = cmfem_parse_double_option(argc, argv, &i, "-a", "--alpha",
                                         &alpha);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-k", "--kappa",
                                         &kappa);
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

   if (kappa < 0.0)
   {
      kappa = (double)(order + 1) * (double)(order + 1);
   }

   // 2. Read the mesh from the given mesh file.
   {
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      CMFEM_DgFeCollection *fec;
      CMFEM_FiniteElementSpace *fespace;
      CMFEM_GridFunction *x;
      CMFEM_VectorFunctionCoefficient *init_x;
      CMFEM_Vector *lambda_values;
      CMFEM_Vector *mu_values;
      CMFEM_PWConstCoefficient *lambda_c;
      CMFEM_PWConstCoefficient *mu_c;
      CMFEM_LinearForm *b;
      CMFEM_BilinearForm *a;
      CMFEM_GSSmoother *M;
      int dim;

      _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list =
         CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_SparseMatrix A = CMFEM_SparseMatrix_Construct();
      _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();

      dim = CMFEM_Mesh_Dimension(mesh);
      if (CMFEM_Mesh_AttributesMax(mesh) < 2 ||
          CMFEM_Mesh_BoundaryAttributesMax(mesh) < 2)
      {
         fprintf(stderr,
                 "\nInput mesh should have at least two materials and two boundary attributes!\n\n");
         CMFEM_ArrayInt_Destroy(&ess_tdof_list);
         CMFEM_SparseMatrix_Destroy(&A);
         CMFEM_Vector_Destroy(&B);
         CMFEM_Vector_Destroy(&X);
         CMFEM_Mesh_Delete(mesh);
         return 3;
      }

      _Alignas(max_align_t) CMFEM_ArrayInt dir_bdr =
         CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));

      // 3. Refine the mesh and convert NURBS geometry if necessary.
      if (ref_levels < 0)
      {
         ref_levels = cmfem_uniform_refinement_levels(5000.0,
                                                      CMFEM_Mesh_GetNE(mesh),
                                                      dim);
      }
      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }
      if (CMFEM_Mesh_HasNURBSext(mesh))
      {
         CMFEM_Mesh_SetCurvature(mesh, order);
      }

      // 4. Define the DG vector finite element space on the mesh.
      fec = CMFEM_DgFeCollection_NewOrderDimBasis(order,
                                                  dim,
                                                  CMFEM_BASIS_GAUSS_LOBATTO);
      fespace = CMFEM_FiniteElementSpace_NewMeshDgVDim(mesh, fec, dim);
      printf("Number of finite element unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(fespace));
      printf("Assembling: ");
      fflush(stdout);

      // 5. Mark attributes 1 and 2 as weak Dirichlet boundaries.
      CMFEM_ArrayInt_Assign(&dir_bdr, 0);
      CMFEM_ArrayInt_Set(&dir_bdr, 0, 1);
      CMFEM_ArrayInt_Set(&dir_bdr, 1, 1);

      // 6. Initialize the DG solution vector with the prescribed displacement.
      x = CMFEM_GridFunction_New(fespace);
      init_x = CMFEM_VectorFunctionCoefficient_New(dim, init_displacement, NULL);
      CMFEM_GridFunction_ProjectCoefficientVfc(x, init_x);

      // 7. Set up the piecewise-constant Lamé coefficients.
      lambda_values = CMFEM_Vector_NewSize(CMFEM_Mesh_AttributesMax(mesh));
      mu_values = CMFEM_Vector_NewSize(CMFEM_Mesh_AttributesMax(mesh));
      CMFEM_Vector_Assign(lambda_values, 1.0);
      CMFEM_Vector_Assign(mu_values, 1.0);
      CMFEM_Vector_Set(lambda_values, 0, 50.0);
      CMFEM_Vector_Set(mu_values, 0, 50.0);
      lambda_c = CMFEM_PWConstCoefficient_New(lambda_values);
      mu_c = CMFEM_PWConstCoefficient_New(mu_values);
      CMFEM_Vector_Delete(lambda_values);
      CMFEM_Vector_Delete(mu_values);

      // 8. Set up the linear form that weakly imposes the Dirichlet data.
      b = CMFEM_LinearForm_New(fespace);
      printf("r.h.s. ... ");
      fflush(stdout);
      CMFEM_LinearForm_AddBdrFaceIntegratorDgeliVfcPwcPwcAi(
         b, init_x, lambda_c, mu_c, alpha, kappa, &dir_bdr);
      CMFEM_LinearForm_Assemble(b);

      // 9. Set up the DG elasticity bilinear form.
      a = CMFEM_BilinearForm_New(fespace);
      CMFEM_BilinearForm_AddDomainIntegratorEi(a, lambda_c, mu_c);
      CMFEM_BilinearForm_AddInteriorFaceIntegratorDgeiPwcPwc(
         a, lambda_c, mu_c, alpha, kappa);
      CMFEM_BilinearForm_AddBdrFaceIntegratorDgeiPwcPwcAi(
         a, lambda_c, mu_c, alpha, kappa, &dir_bdr);

      // 10. Assemble the bilinear form and the corresponding linear system.
      printf("matrix ... ");
      fflush(stdout);
      CMFEM_BilinearForm_Assemble(a);
      CMFEM_BilinearForm_FormLinearSystemSm(a, &ess_tdof_list, x, b, &A, &X, &B);
      printf("done.\n");

      // 11. Solve the symmetric or non-symmetric system with a GS
      //     preconditioner.
      M = CMFEM_GSSmoother_NewSm(&A);
      if (alpha == -1.0)
      {
         CMFEM_PCGSmGs(&A, M, &B, &X, 3, 5000, 1e-12, 0.0);
      }
      else
      {
         CMFEM_GMRESSmGs(&A, M, &B, &X, 3, 5000, 100, 1e-12, 0.0);
      }
      CMFEM_GSSmoother_Delete(M);

      // 12. Recover the displacement field and save the deformed mesh.
      CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
      CMFEM_Mesh_SetNodalFESpace(mesh, fespace);
      CMFEM_Mesh_AddDisplacementToNodes(mesh, x);
      CMFEM_GridFunction_Scale(x, -1.0);
      CMFEM_Mesh_Print(mesh, "displaced.mesh", 8);
      CMFEM_GridFunction_Save(x, "sol.gf", 8);

      // 13. Send the displaced configuration to GLVis.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
      }

      CMFEM_BilinearForm_Delete(a);
      CMFEM_LinearForm_Delete(b);
      CMFEM_PWConstCoefficient_Delete(mu_c);
      CMFEM_PWConstCoefficient_Delete(lambda_c);
      CMFEM_VectorFunctionCoefficient_Delete(init_x);
      CMFEM_GridFunction_Delete(x);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_DgFeCollection_Delete(fec);
      CMFEM_ArrayInt_Destroy(&dir_bdr);
      CMFEM_ArrayInt_Destroy(&ess_tdof_list);
      CMFEM_SparseMatrix_Destroy(&A);
      CMFEM_Vector_Destroy(&B);
      CMFEM_Vector_Destroy(&X);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
