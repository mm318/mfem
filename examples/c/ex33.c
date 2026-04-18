//                              MFEM C Example 33
//
// Description: This example solves a fractional Laplace problem by first
//              reducing the integer-order part and then summing a rational
//              partial-fraction expansion of shifted Helmholtz solves.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

static const double cmfem_pi_ = 3.14159265358979323846;

static double fractional_source(const CMFEM_Vector *x, void *context)
{
   const double alpha = *(const double *)context;
   const int dim = CMFEM_Vector_Size(x);
   double value = 1.0;
   int i;

   for (i = 0; i < dim; i++)
   {
      value *= sin(cmfem_pi_ * CMFEM_Vector_Get(x, i));
   }

   return pow(dim * cmfem_pi_ * cmfem_pi_, alpha) * value;
}

static double manufactured_solution(const CMFEM_Vector *x, void *context)
{
   const int dim = CMFEM_Vector_Size(x);
   double value = 1.0;
   int i;

   (void)context;
   for (i = 0; i < dim; i++)
   {
      value *= sin(cmfem_pi_ * CMFEM_Vector_Get(x, i));
   }
   return value;
}

static void print_verification_summary(const char *mesh_file,
                                       int dim,
                                       double l2_error)
{
   const char *manufactured_solution_text = "sin(pi x)";
   const char *expected_mesh = "inline_segment.mesh";

   if (dim == 2)
   {
      manufactured_solution_text = "sin(pi x) sin(pi y)";
      expected_mesh = "inline_quad.mesh";
   }
   else if (dim >= 3)
   {
      manufactured_solution_text = "sin(pi x) sin(pi y) sin(pi z)";
      expected_mesh = "inline_hex.mesh";
   }

   printf("\n");
   printf("================================================================================\n\n");
   printf("Solution Verification in %dD\n\n", dim);
   printf("Manufactured solution : %s\n", manufactured_solution_text);
   printf("Expected mesh         : %s\n", expected_mesh);
   printf("Your mesh             : %s\n", mesh_file);
   printf("L2 error              : %.8g\n\n", l2_error);
   printf("================================================================================\n");
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int num_refs = 3;
   double alpha = 0.5;
   int visualization = 1;
   int verification = 0;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-r", "--refs",
                                      &num_refs);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-alpha", "--alpha",
                                         &alpha);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-ver", "--verification",
                                       "-no-ver", "--no-verification",
                                       &verification);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Compute the rational partial-fraction data for the non-integer part.
   {
      const int power_of_laplace = (int)floor(alpha);
      double exponent_to_approximate = alpha - power_of_laplace;
      int integer_order = fabs(exponent_to_approximate) <= 1.0e-12;
      int coeff_count = 0;
      double *coeffs = NULL;
      double *poles = NULL;
      CMFEM_Mesh *mesh = NULL;
      CMFEM_H1FeCollection *fec = NULL;
      CMFEM_FiniteElementSpace *fespace = NULL;
      CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
      CMFEM_GridFunction *u = NULL;
      CMFEM_GridFunction *x = NULL;
      CMFEM_GridFunction *g = NULL;
      CMFEM_LinearForm *b = NULL;
      CMFEM_ConstantCoefficient *one = NULL;
      CMFEM_FunctionCoefficient *f = NULL;
      int dim;

      if (!integer_order)
      {
         printf("Approximating the fractional exponent %.12g\n",
                exponent_to_approximate);
         CMFEM_ComputePartialFractionApproximation(&exponent_to_approximate,
                                                   &coeff_count,
                                                   &coeffs,
                                                   &poles);
         alpha = exponent_to_approximate + power_of_laplace;
      }
      else
      {
         printf("Treating integer order PDE.\n");
      }

      // 3. Read and refine the mesh.
      mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      dim = CMFEM_Mesh_Dimension(mesh);
      for (i = 0; i < num_refs; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      // 4. Define the scalar H1 finite element space and essential boundary
      //    conditions.
      fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, fec);
      printf("Number of degrees of freedom: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

      if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
      {
         _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
            CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
         CMFEM_ArrayInt_Assign(&ess_bdr, 1);
         CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace,
                                                       &ess_bdr,
                                                       &ess_tdof_list);
         CMFEM_ArrayInt_Destroy(&ess_bdr);
      }

      // 5. Define the initial right-hand side and solution grid functions.
      one = CMFEM_ConstantCoefficient_New(1.0);
      if (verification)
      {
         f = CMFEM_FunctionCoefficient_New(fractional_source, &alpha);
         b = CMFEM_LinearForm_New(fespace);
         CMFEM_LinearForm_AddDomainIntegratorDliFc(b, f);
      }
      else
      {
         b = CMFEM_LinearForm_New(fespace);
         CMFEM_LinearForm_AddDomainIntegratorDliCc(b, one);
      }
      CMFEM_LinearForm_Assemble(b);

      u = CMFEM_GridFunction_New(fespace);
      x = CMFEM_GridFunction_New(fespace);
      g = CMFEM_GridFunction_New(fespace);
      CMFEM_GridFunction_Assign(u, 0.0);
      CMFEM_GridFunction_Assign(x, 0.0);
      CMFEM_GridFunction_Assign(g, 0.0);

      // 6. Solve the integer-order reduction (-Delta)^N g = f.
      if (power_of_laplace > 0)
      {
         CMFEM_BilinearForm *k = CMFEM_BilinearForm_New(fespace);
         CMFEM_BilinearForm *m = CMFEM_BilinearForm_New(fespace);
         CMFEM_GSSmoother *smoother = NULL;
         _Alignas(max_align_t) CMFEM_SparseMatrix A =
            CMFEM_SparseMatrix_Construct();
         _Alignas(max_align_t) CMFEM_SparseMatrix mass =
            CMFEM_SparseMatrix_Construct();
         _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
         _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();

         CMFEM_BilinearForm_AddDomainIntegratorDiCc(k, one);
         CMFEM_BilinearForm_Assemble(k);
         CMFEM_BilinearForm_AddDomainIntegratorMiCc(m, one);
         CMFEM_BilinearForm_Assemble(m);
         CMFEM_BilinearForm_FormSystemMatrixSm(m, &ess_tdof_list, &mass);
         CMFEM_BilinearForm_FormLinearSystemSm(k,
                                               &ess_tdof_list,
                                               g,
                                               b,
                                               &A,
                                               &X,
                                               &B);
         smoother = CMFEM_GSSmoother_NewSm(&A);

         for (i = 0; i < power_of_laplace; i++)
         {
            CMFEM_PCGSmGs(&A, smoother, &B, &X, 3, 300, 1.0e-12, 0.0);
            CMFEM_BilinearForm_RecoverFEMSolution(k, &X, b, g);
            CMFEM_SparseMatrix_Mult(&mass, &X, &B);

            if (integer_order && i == power_of_laplace - 1)
            {
               CMFEM_GridFunction_Add(u, g);
            }
         }

         if (!integer_order)
         {
            const CMFEM_SparseMatrix *restriction =
               CMFEM_FiniteElementSpace_GetRestrictionMatrixSm(fespace);
            _Alignas(max_align_t) CMFEM_Vector rhs_dofs =
               CMFEM_Vector_ConstructSize(CMFEM_FiniteElementSpace_GetVSize(
                                             fespace));

            if (restriction != NULL)
            {
               CMFEM_SparseMatrix_MultTranspose(restriction, &B, &rhs_dofs);
               CMFEM_LinearForm_SetFromVec(b, &rhs_dofs);
            }
            else
            {
               CMFEM_LinearForm_SetFromVec(b, &B);
            }

            CMFEM_Vector_Destroy(&rhs_dofs);
         }

         CMFEM_GSSmoother_Delete(smoother);
         CMFEM_Vector_Destroy(&B);
         CMFEM_Vector_Destroy(&X);
         CMFEM_SparseMatrix_Destroy(&mass);
         CMFEM_SparseMatrix_Destroy(&A);
         CMFEM_BilinearForm_Delete(m);
         CMFEM_BilinearForm_Delete(k);
      }

      // 7. Solve the fractional PDE as a sum of shifted Helmholtz solves.
      if (!integer_order)
      {
         for (i = 0; i < coeff_count; i++)
         {
            CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
            CMFEM_ConstantCoefficient *shift =
               CMFEM_ConstantCoefficient_New(-poles[i]);
            CMFEM_GSSmoother *smoother = NULL;
            _Alignas(max_align_t) CMFEM_SparseMatrix A =
               CMFEM_SparseMatrix_Construct();
            _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
            _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();

            printf("\nSolving PDE -Delta u + %.12g u = %.12g g\n",
                   -poles[i], coeffs[i]);

            CMFEM_GridFunction_Assign(x, 0.0);
            CMFEM_BilinearForm_AddDomainIntegratorDiCc(a, one);
            CMFEM_BilinearForm_AddDomainIntegratorMiCc(a, shift);
            CMFEM_BilinearForm_Assemble(a);
            CMFEM_BilinearForm_FormLinearSystemSm(a,
                                                  &ess_tdof_list,
                                                  x,
                                                  b,
                                                  &A,
                                                  &X,
                                                  &B);

            smoother = CMFEM_GSSmoother_NewSm(&A);
            CMFEM_PCGSmGs(&A, smoother, &B, &X, 3, 300, 1.0e-12, 0.0);
            CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
            CMFEM_GridFunction_Scale(x, coeffs[i]);
            CMFEM_GridFunction_Add(u, x);

            CMFEM_GSSmoother_Delete(smoother);
            CMFEM_Vector_Destroy(&B);
            CMFEM_Vector_Destroy(&X);
            CMFEM_SparseMatrix_Destroy(&A);
            CMFEM_ConstantCoefficient_Delete(shift);
            CMFEM_BilinearForm_Delete(a);
         }
      }

      // 8. Optionally verify the final solution.
      if (verification)
      {
         CMFEM_FunctionCoefficient *solution =
            CMFEM_FunctionCoefficient_New(manufactured_solution, NULL);
         const double l2_error =
            CMFEM_GridFunction_ComputeL2ErrorFc(u, solution);
         print_verification_summary(mesh_file, dim, l2_error);
         CMFEM_FunctionCoefficient_Delete(solution);
      }

      // 9. Send the final solution to GLVis.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, u, "localhost", 19916);
      }

      CMFEM_FreeDoubles(poles);
      CMFEM_FreeDoubles(coeffs);
      CMFEM_GridFunction_Delete(g);
      CMFEM_GridFunction_Delete(x);
      CMFEM_GridFunction_Delete(u);
      CMFEM_LinearForm_Delete(b);
      if (f != NULL) { CMFEM_FunctionCoefficient_Delete(f); }
      CMFEM_ConstantCoefficient_Delete(one);
      CMFEM_ArrayInt_Destroy(&ess_tdof_list);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_H1FeCollection_Delete(fec);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
