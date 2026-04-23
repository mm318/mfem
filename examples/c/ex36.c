//                              MFEM C Example 36
//
// Description: This example solves a bound-constrained obstacle problem using
//              a proximal Galerkin iteration on coupled H1/L2 spaces.

#include "common.h"
#include "cmfem.h"
#include "adapters/ex36/obstacle_assembly.h"

#include <math.h>
#include <stdio.h>

static double zero_function(const CMFEM_Vector *x, void *context)
{
   (void)x;
   (void)context;
   return 0.0;
}

static void copy_vector(CMFEM_Vector *dst, const CMFEM_Vector *src)
{
   int i;
   for (i = 0; i < CMFEM_Vector_Size(src); i++)
   {
      CMFEM_Vector_Set(dst, i, CMFEM_Vector_Get(src, i));
   }
}

static void pack_block(CMFEM_Vector *dst, int offset, const CMFEM_Vector *src)
{
   int i;
   for (i = 0; i < CMFEM_Vector_Size(src); i++)
   {
      CMFEM_Vector_Set(dst, offset + i, CMFEM_Vector_Get(src, i));
   }
}

static void unpack_block(CMFEM_Vector *dst, const CMFEM_Vector *src, int offset)
{
   int i;
   for (i = 0; i < CMFEM_Vector_Size(dst); i++)
   {
      CMFEM_Vector_Set(dst, i, CMFEM_Vector_Get(src, offset + i));
   }
}

static void scale_mesh_point(const CMFEM_Vector *x, CMFEM_Vector *value,
                             void *context)
{
   const double scale = *(const double *)context;
   int i;
   for (i = 0; i < CMFEM_Vector_Size(x); i++)
   {
      CMFEM_Vector_Set(value, i, CMFEM_Vector_Get(x, i) / scale);
   }
}

static double spherical_obstacle(const CMFEM_Vector *pt, void *context)
{
   const double x = CMFEM_Vector_Get(pt, 0);
   const double y = CMFEM_Vector_Get(pt, 1);
   const double r = sqrt(x * x + y * y);
   const double r0 = 0.5;
   const double beta = 0.9;
   const double b = r0 * beta;
   const double tmp = sqrt(r0 * r0 - b * b);
   const double B = tmp + b * b / tmp;
   const double C = -b / tmp;

   (void)context;
   if (r > b)
   {
      return B + r * C;
   }
   return sqrt(r0 * r0 - r * r);
}

static double exact_solution_obstacle(const CMFEM_Vector *pt, void *context)
{
   const double x = CMFEM_Vector_Get(pt, 0);
   const double y = CMFEM_Vector_Get(pt, 1);
   const double r = sqrt(x * x + y * y);
   const double r0 = 0.5;
   const double a = 0.348982574111686;
   const double A = -0.340129705945858;

   (void)context;
   if (r > a)
   {
      return A * log(r);
   }
   return sqrt(r0 * r0 - r * r);
}

static void exact_solution_gradient_obstacle(const CMFEM_Vector *pt,
                                             CMFEM_Vector *grad,
                                             void *context)
{
   const double x = CMFEM_Vector_Get(pt, 0);
   const double y = CMFEM_Vector_Get(pt, 1);
   const double r = sqrt(x * x + y * y);
   const double r0 = 0.5;
   const double a = 0.348982574111686;
   const double A = -0.340129705945858;

   (void)context;
   if (r > a)
   {
      CMFEM_Vector_Set(grad, 0, A * x / (r * r));
      CMFEM_Vector_Set(grad, 1, A * y / (r * r));
   }
   else
   {
      CMFEM_Vector_Set(grad, 0, -x / sqrt(r0 * r0 - r * r));
      CMFEM_Vector_Set(grad, 1, -y / sqrt(r0 * r0 - r * r));
   }
}

static double initial_condition(const CMFEM_Vector *x, void *context)
{
   double rr = 0.0;
   int i;

   (void)context;
   for (i = 0; i < CMFEM_Vector_Size(x); i++)
   {
      const double xi = CMFEM_Vector_Get(x, i);
      rr += xi * xi;
   }
   return 1.0 - rr;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int order = 1;
   int max_it = 10;
   int ref_levels = 3;
   double alpha = 1.0;
   double tol = 1e-5;
   int visualization = 1;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-r", "--refs",
                                      &ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-mi", "--max-it",
                                      &max_it);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-tol", "--tol",
                                         &tol);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-step", "--step",
                                         &alpha);
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

   // 2. Read and postprocess the curved unit-disc mesh.
   {
      const char *mesh_file = CMFEM_ExamplesDataPath("disc-nurbs.mesh");
      const double mesh_scale = 2.0 * sqrt(2.0);
      const double log_min_value = -36.0;
      const double exp_min_value = 0.0;
      const double exp_max_value = 1.0e6;
      const double spectrum_shift = -1.0e-6;
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      int dim = CMFEM_Mesh_Dimension(mesh);
      int curvature_order;
      CMFEM_H1FeCollection *h1_fec;
      CMFEM_L2FeCollection *l2_fec;
      CMFEM_FiniteElementSpace *h1_fes;
      CMFEM_FiniteElementSpace *l2_fes;
      int h1_size;
      int l2_size;
      int h1_true_size;
      int l2_true_size;
      int total_iterations = 0;
      double increment_u = 0.1;
      int outer_iterations = 0;
      _Alignas(max_align_t) CMFEM_ArrayInt offsets =
         CMFEM_ArrayInt_ConstructSize(3);
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
         CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
      _Alignas(max_align_t) CMFEM_Vector x = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector rhs = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector rhs_h1 = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector rhs_h1_tmp = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector rhs_l2 = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector u_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector delta_psi_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector u_old_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector psi_old_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector psi_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector u_prev_newton_vec =
         CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector tmp_u_diff = CMFEM_Vector_Construct();
      CMFEM_GridFunction *u_gf;
      CMFEM_GridFunction *delta_psi_gf;
      CMFEM_GridFunction *u_old_gf;
      CMFEM_GridFunction *psi_old_gf;
      CMFEM_GridFunction *psi_gf;
      CMFEM_GridFunction *u_tmp_gf;
      CMFEM_GridFunction *u_alt_gf;
      CMFEM_FunctionCoefficient *zero_cf;
      CMFEM_FunctionCoefficient *exact_cf;
      CMFEM_VectorFunctionCoefficient *exact_grad_cf;
      CMFEM_FunctionCoefficient *initial_cf;
      CMFEM_FunctionCoefficient *obstacle_cf;

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }
      curvature_order = order > 2 ? order : 2;
      CMFEM_Mesh_SetCurvature(mesh, curvature_order);
      CMFEM_Mesh_Transform(mesh, scale_mesh_point, (void *)&mesh_scale);

      // 3. Define the coupled H1/L2 finite element spaces.
      h1_fec = CMFEM_H1FeCollection_NewOrderDim(order + 1, dim);
      l2_fec = CMFEM_L2FeCollection_NewOrderDim(order - 1, dim);
      h1_fes = CMFEM_FiniteElementSpace_NewMeshH1(mesh, h1_fec);
      l2_fes = CMFEM_FiniteElementSpace_NewMeshL2Vdim(mesh, l2_fec, 1);

      h1_size = CMFEM_FiniteElementSpace_GetVSize(h1_fes);
      l2_size = CMFEM_FiniteElementSpace_GetVSize(l2_fes);
      h1_true_size = CMFEM_FiniteElementSpace_GetTrueVSize(h1_fes);
      l2_true_size = CMFEM_FiniteElementSpace_GetTrueVSize(l2_fes);

      printf("Number of H1 finite element unknowns: %d\n", h1_true_size);
      printf("Number of L2 finite element unknowns: %d\n", l2_true_size);

      CMFEM_ArrayInt_Set(&offsets, 0, 0);
      CMFEM_ArrayInt_Set(&offsets, 1, h1_size);
      CMFEM_ArrayInt_Set(&offsets, 2, h1_size + l2_size);
      CMFEM_ArrayInt_Assign(&ess_bdr, 1);

      x = CMFEM_Vector_ConstructSize(h1_size + l2_size);
      rhs = CMFEM_Vector_ConstructSize(h1_size + l2_size);
      rhs_h1 = CMFEM_Vector_ConstructSize(h1_size);
      rhs_h1_tmp = CMFEM_Vector_ConstructSize(h1_size);
      rhs_l2 = CMFEM_Vector_ConstructSize(l2_size);
      u_vec = CMFEM_Vector_ConstructSize(h1_size);
      delta_psi_vec = CMFEM_Vector_ConstructSize(l2_size);
      u_old_vec = CMFEM_Vector_ConstructSize(h1_size);
      psi_old_vec = CMFEM_Vector_ConstructSize(l2_size);
      psi_vec = CMFEM_Vector_ConstructSize(l2_size);
      u_prev_newton_vec = CMFEM_Vector_ConstructSize(h1_size);
      tmp_u_diff = CMFEM_Vector_ConstructSize(h1_size);

      CMFEM_Vector_Assign(&x, 0.0);
      CMFEM_Vector_Assign(&rhs, 0.0);
      CMFEM_Vector_Assign(&rhs_h1, 0.0);
      CMFEM_Vector_Assign(&rhs_h1_tmp, 0.0);
      CMFEM_Vector_Assign(&rhs_l2, 0.0);
      CMFEM_Vector_Assign(&u_vec, 0.0);
      CMFEM_Vector_Assign(&delta_psi_vec, 0.0);
      CMFEM_Vector_Assign(&u_old_vec, 0.0);
      CMFEM_Vector_Assign(&psi_old_vec, 0.0);
      CMFEM_Vector_Assign(&psi_vec, 0.0);
      CMFEM_Vector_Assign(&u_prev_newton_vec, 0.0);
      CMFEM_Vector_Assign(&tmp_u_diff, 0.0);

      // 4. Define the solution grid functions and scalar coefficients.
      u_gf = CMFEM_GridFunction_New(h1_fes);
      delta_psi_gf = CMFEM_GridFunction_New(l2_fes);
      u_old_gf = CMFEM_GridFunction_New(h1_fes);
      psi_old_gf = CMFEM_GridFunction_New(l2_fes);
      psi_gf = CMFEM_GridFunction_New(l2_fes);
      u_tmp_gf = CMFEM_GridFunction_New(h1_fes);
      u_alt_gf = CMFEM_GridFunction_New(l2_fes);

      CMFEM_GridFunction_Assign(u_gf, 0.0);
      CMFEM_GridFunction_Assign(delta_psi_gf, 0.0);
      CMFEM_GridFunction_Assign(u_old_gf, 0.0);
      CMFEM_GridFunction_Assign(psi_old_gf, 0.0);
      CMFEM_GridFunction_Assign(psi_gf, 0.0);
      CMFEM_GridFunction_Assign(u_tmp_gf, 0.0);
      CMFEM_GridFunction_Assign(u_alt_gf, 0.0);

      zero_cf = CMFEM_FunctionCoefficient_New(zero_function, NULL);
      exact_cf = CMFEM_FunctionCoefficient_New(exact_solution_obstacle, NULL);
      exact_grad_cf = CMFEM_VectorFunctionCoefficient_New(
                         dim,
                         exact_solution_gradient_obstacle,
                         NULL);
      initial_cf = CMFEM_FunctionCoefficient_New(initial_condition, NULL);
      obstacle_cf = CMFEM_FunctionCoefficient_New(spherical_obstacle, NULL);

      // 5. Initialize the primal solution and slack variable psi_h.
      CMFEM_GridFunction_ProjectCoefficientFc(u_gf, initial_cf);
      CMFEM_GridFunction_CopyToVec(u_gf, &u_vec);
      copy_vector(&u_old_vec, &u_vec);
      pack_block(&x, 0, &u_vec);
      CMFEM_GridFunction_SetFromVec(u_old_gf, &u_old_vec);

      CMFEM_ProjectObstacleLogGfFc(psi_gf,
                                   u_gf,
                                   obstacle_cf,
                                   log_min_value);
      CMFEM_GridFunction_CopyToVec(psi_gf, &psi_vec);
      copy_vector(&psi_old_vec, &psi_vec);
      CMFEM_GridFunction_SetFromVec(psi_old_gf, &psi_old_vec);

      // 6. Run the proximal Galerkin outer loop and the inner Newton solves.
      for (i = 0; i < max_it; i++)
      {
         int j;
         outer_iterations = i + 1;
         copy_vector(&u_prev_newton_vec, &u_old_vec);

         printf("\nOUTER ITERATION %d\n", i + 1);

         for (j = 0; j < 10; j++)
         {
            CMFEM_GridFunctionCoefficient *psi_old_cf;
            CMFEM_GridFunctionCoefficient *psi_cf;
            CMFEM_LinearForm *b0_lf;
            CMFEM_LinearForm *b0_tmp_lf;
            CMFEM_ConstantCoefficient *alpha_cf;
            CMFEM_BilinearForm *a00;
            CMFEM_MixedBilinearForm *a10;
            CMFEM_SparseMatrix *a00_matrix;
            CMFEM_SparseMatrix *a10_matrix;
            CMFEM_SparseMatrix *a01_matrix;
            CMFEM_SparseMatrix *a11_matrix;
            CMFEM_BlockOperator *a_block;
            CMFEM_BlockDiagonalPreconditioner *prec;
            CMFEM_GSSmoother *a00_prec;
            CMFEM_GSSmoother *a11_prec;
            double newton_update_size;

            total_iterations++;

            psi_old_cf = CMFEM_GridFunctionCoefficient_NewGf(psi_old_gf);
            psi_cf = CMFEM_GridFunctionCoefficient_NewGf(psi_gf);

            b0_lf = CMFEM_LinearForm_New(h1_fes);
            CMFEM_LinearForm_AddDomainIntegratorDliGfc(b0_lf, psi_old_cf);
            CMFEM_LinearForm_Assemble(b0_lf);
            CMFEM_LinearForm_CopyToVector(b0_lf, &rhs_h1);

            b0_tmp_lf = CMFEM_LinearForm_New(h1_fes);
            CMFEM_LinearForm_AddDomainIntegratorDliGfc(b0_tmp_lf, psi_cf);
            CMFEM_LinearForm_Assemble(b0_tmp_lf);
            CMFEM_LinearForm_CopyToVector(b0_tmp_lf, &rhs_h1_tmp);
            CMFEM_Vector_AddScaled(&rhs_h1, -1.0, &rhs_h1_tmp);

            CMFEM_LinearForm_Delete(b0_tmp_lf);
            CMFEM_LinearForm_Delete(b0_lf);
            CMFEM_GridFunctionCoefficient_Delete(psi_cf);
            CMFEM_GridFunctionCoefficient_Delete(psi_old_cf);

            CMFEM_AssembleObstacleLfGfFc(l2_fes,
                                         psi_gf,
                                         obstacle_cf,
                                         &rhs_l2,
                                         exp_min_value,
                                         exp_max_value);

            alpha_cf = CMFEM_ConstantCoefficient_New(alpha);
            a00 = CMFEM_BilinearForm_New(h1_fes);
            CMFEM_BilinearForm_SetDiagonalPolicyOne(a00);
            CMFEM_BilinearForm_AddDomainIntegratorDiCc(a00, alpha_cf);
            CMFEM_BilinearForm_Assemble(a00);
            CMFEM_BilinearForm_EliminateEssentialBCAiVecVec(a00,
                                                            &ess_bdr,
                                                            &u_vec,
                                                            &rhs_h1);
            CMFEM_BilinearForm_Finalize(a00);
            a00_matrix = CMFEM_BilinearForm_SpMat(a00);

            a10 = CMFEM_MixedBilinearForm_New(h1_fes, l2_fes);
            CMFEM_MixedBilinearForm_AddDomainIntegratorMsi(a10);
            CMFEM_MixedBilinearForm_Assemble(a10);
            CMFEM_MixedBilinearForm_EliminateTrialEssentialBCAiVecVec(a10,
                                                                      &ess_bdr,
                                                                      &u_vec,
                                                                      &rhs_l2);
            CMFEM_MixedBilinearForm_Finalize(a10);
            a10_matrix = CMFEM_MixedBilinearForm_SpMat(a10);
            a01_matrix = CMFEM_TransposeSm(a10_matrix);
            a11_matrix = CMFEM_AssembleObstacleHessianSmGf(l2_fes,
                                                           psi_gf,
                                                           order,
                                                           spectrum_shift,
                                                           exp_min_value,
                                                           exp_max_value);

            pack_block(&rhs, 0, &rhs_h1);
            pack_block(&rhs, h1_size, &rhs_l2);

            a_block = CMFEM_BlockOperator_NewAiAi(&offsets, &offsets);
            CMFEM_BlockOperator_SetBlockSm(a_block, 0, 0, a00_matrix);
            CMFEM_BlockOperator_SetBlockSm(a_block, 1, 0, a10_matrix);
            CMFEM_BlockOperator_SetBlockSm(a_block, 0, 1, a01_matrix);
            CMFEM_BlockOperator_SetBlockSm(a_block, 1, 1, a11_matrix);

            a00_prec = CMFEM_GSSmoother_NewSm(a00_matrix);
            a11_prec = CMFEM_GSSmoother_NewSm(a11_matrix);
            prec = CMFEM_BlockDiagonalPreconditioner_NewAi(&offsets);
            CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockGs(prec, 0, a00_prec);
            CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockGs(prec, 1, a11_prec);

            CMFEM_GMRESBopBdp(a_block, prec, &rhs, &x, 0, 10000, 500, 1e-12, 0.0);

            unpack_block(&u_vec, &x, 0);
            unpack_block(&delta_psi_vec, &x, h1_size);
            CMFEM_GridFunction_SetFromVec(u_gf, &u_vec);
            CMFEM_GridFunction_SetFromVec(delta_psi_gf, &delta_psi_vec);

            copy_vector(&tmp_u_diff, &u_prev_newton_vec);
            CMFEM_Vector_AddScaled(&tmp_u_diff, -1.0, &u_vec);
            CMFEM_GridFunction_SetFromVec(u_tmp_gf, &tmp_u_diff);
            newton_update_size =
               CMFEM_GridFunction_ComputeL2ErrorFc(u_tmp_gf, zero_cf);
            copy_vector(&u_prev_newton_vec, &u_vec);

            CMFEM_Vector_AddScaled(&psi_vec, 1.0, &delta_psi_vec);
            CMFEM_GridFunction_SetFromVec(psi_gf, &psi_vec);

            if (visualization)
            {
               CMFEM_SendSolutionToGLVis(mesh, u_gf, "localhost", 19916);
            }

            CMFEM_BlockDiagonalPreconditioner_Delete(prec);
            CMFEM_GSSmoother_Delete(a11_prec);
            CMFEM_GSSmoother_Delete(a00_prec);
            CMFEM_BlockOperator_Delete(a_block);
            CMFEM_SparseMatrix_Delete(a11_matrix);
            CMFEM_SparseMatrix_Delete(a01_matrix);
            CMFEM_MixedBilinearForm_Delete(a10);
            CMFEM_BilinearForm_Delete(a00);
            CMFEM_ConstantCoefficient_Delete(alpha_cf);

            if (newton_update_size < increment_u)
            {
               break;
            }
         }

         copy_vector(&tmp_u_diff, &u_vec);
         CMFEM_Vector_AddScaled(&tmp_u_diff, -1.0, &u_old_vec);
         CMFEM_GridFunction_SetFromVec(u_tmp_gf, &tmp_u_diff);
         increment_u = CMFEM_GridFunction_ComputeL2ErrorFc(u_tmp_gf, zero_cf);

         printf("Number of Newton iterations = %d\n", j + 1);
         printf("Increment (|| u_h - u_h_prvs||) = %.8g\n", increment_u);

         copy_vector(&u_old_vec, &u_vec);
         copy_vector(&psi_old_vec, &psi_vec);
         CMFEM_GridFunction_SetFromVec(u_old_gf, &u_old_vec);
         CMFEM_GridFunction_SetFromVec(psi_old_gf, &psi_old_vec);

         if (increment_u < tol || i == max_it - 1)
         {
            break;
         }

         printf("H1-error  (|| u - u_h^k||)       = %.8g\n",
                CMFEM_GridFunction_ComputeH1ErrorFcVfc(u_gf,
                                                       exact_cf,
                                                       exact_grad_cf));
      }

      printf("\n Outer iterations: %d\n", outer_iterations);
      printf(" Total iterations: %d\n", total_iterations);
      printf(" Total dofs:       %d\n", h1_true_size + l2_true_size);

      // 7. Report the final discretization errors.
      CMFEM_ProjectObstacleExpGfFc(u_alt_gf,
                                   psi_gf,
                                   obstacle_cf,
                                   exp_min_value,
                                   exp_max_value);
      printf("\n Final L2-error (|| u - u_h||)          = %.8g\n",
             CMFEM_GridFunction_ComputeL2ErrorFc(u_gf, exact_cf));
      printf(" Final H1-error (|| u - u_h||)          = %.8g\n",
             CMFEM_GridFunction_ComputeH1ErrorFcVfc(u_gf,
                                                    exact_cf,
                                                    exact_grad_cf));
      printf(" Final L2-error (|| u - phi - exp(psi_h)||) = %.8g\n",
             CMFEM_GridFunction_ComputeL2ErrorFc(u_alt_gf, exact_cf));

      CMFEM_FunctionCoefficient_Delete(obstacle_cf);
      CMFEM_FunctionCoefficient_Delete(initial_cf);
      CMFEM_VectorFunctionCoefficient_Delete(exact_grad_cf);
      CMFEM_FunctionCoefficient_Delete(exact_cf);
      CMFEM_FunctionCoefficient_Delete(zero_cf);
      CMFEM_GridFunction_Delete(u_alt_gf);
      CMFEM_GridFunction_Delete(u_tmp_gf);
      CMFEM_GridFunction_Delete(psi_gf);
      CMFEM_GridFunction_Delete(psi_old_gf);
      CMFEM_GridFunction_Delete(u_old_gf);
      CMFEM_GridFunction_Delete(delta_psi_gf);
      CMFEM_GridFunction_Delete(u_gf);
      CMFEM_Vector_Destroy(&tmp_u_diff);
      CMFEM_Vector_Destroy(&u_prev_newton_vec);
      CMFEM_Vector_Destroy(&psi_vec);
      CMFEM_Vector_Destroy(&psi_old_vec);
      CMFEM_Vector_Destroy(&u_old_vec);
      CMFEM_Vector_Destroy(&delta_psi_vec);
      CMFEM_Vector_Destroy(&u_vec);
      CMFEM_Vector_Destroy(&rhs_l2);
      CMFEM_Vector_Destroy(&rhs_h1_tmp);
      CMFEM_Vector_Destroy(&rhs_h1);
      CMFEM_Vector_Destroy(&rhs);
      CMFEM_Vector_Destroy(&x);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_ArrayInt_Destroy(&offsets);
      CMFEM_FiniteElementSpace_Delete(l2_fes);
      CMFEM_FiniteElementSpace_Delete(h1_fes);
      CMFEM_L2FeCollection_Delete(l2_fec);
      CMFEM_H1FeCollection_Delete(h1_fec);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
