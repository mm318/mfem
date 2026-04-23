//                              MFEM C Example 40
//
// Description: This example solves the eikonal equation by a proximal
//              Galerkin iteration posed as a sequence of mixed RT/L2 problems.

#include "common.h"
#include "cmfem.h"
#include "adapters/ex40/isomorphism_assembly.h"

#include <math.h>
#include <stdio.h>

static double zero_fun(const CMFEM_Vector *x, void *context)
{
   (void)x;
   (void)context;
   return 0.0;
}

static double max_double(double a, double b)
{
   return a > b ? a : b;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int max_it = 5;
   int ref_levels = 3;
   double alpha = 1.0;
   double growth_rate = 1.0;
   double newton_scaling = 0.8;
   double eps = 1e-6;
   double tol = 1e-4;
   int visualization = 1;
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

      parsed = cmfem_parse_double_option(argc, argv, &i, "-gr",
                                         "--growth-rate", &growth_rate);
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

   // 2. Read and refine the mesh.
   {
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      CMFEM_RtFeCollection *rt_fec;
      CMFEM_L2FeCollection *l2_fec;
      CMFEM_FiniteElementSpace *rt_fes;
      CMFEM_FiniteElementSpace *l2_fes;
      CMFEM_MixedBilinearForm *a10;
      CMFEM_SparseMatrix *a10_matrix;
      CMFEM_SparseMatrix *a01_matrix;
      CMFEM_FunctionCoefficient *zero_cf;
      CMFEM_GridFunction *delta_psi_gf;
      CMFEM_GridFunction *u_gf;
      CMFEM_GridFunction *psi_old_gf;
      CMFEM_GridFunction *psi_gf;
      CMFEM_GridFunction *u_old_gf;
      CMFEM_GridFunction *u_tmp_gf;
      int dim;
      int sdim;
      int curvature_order;
      int rt_size;
      int l2_size;
      int outer_iterations = 0;
      int total_iterations = 0;
      double increment_u = 0.1;
      double alpha_k = alpha;

      _Alignas(max_align_t) CMFEM_ArrayInt offsets =
         CMFEM_ArrayInt_ConstructSize(3);
      _Alignas(max_align_t) CMFEM_Vector x = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector rhs = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector delta_psi_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector u_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector psi_old_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector psi_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector u_old_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector u_prev_newton_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector b0_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector b1_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector tmp_l2_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector tmp_l2_vec2 = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector tmp_u_diff = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector a00_diag = CMFEM_Vector_Construct();

      dim = CMFEM_Mesh_Dimension(mesh);
      sdim = CMFEM_Mesh_SpaceDimension(mesh);

      if (CMFEM_Mesh_BoundaryAttributesSize(mesh) == 0)
      {
         fprintf(stderr,
                 "This example does not support meshes without boundary attributes.\n");
         CMFEM_Mesh_Delete(mesh);
         return 1;
      }

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      curvature_order = order > 2 ? order : 2;
      CMFEM_Mesh_SetCurvature(mesh, curvature_order);

      // 3. Define the RT/L2 spaces and solution storage.
      rt_fec = CMFEM_RtFeCollection_NewOrderDim(order, dim);
      l2_fec = CMFEM_L2FeCollection_NewOrderDim(order, dim);
      rt_fes = CMFEM_FiniteElementSpace_NewMeshRt(mesh, rt_fec);
      l2_fes = CMFEM_FiniteElementSpace_NewMeshL2Vdim(mesh, l2_fec, 1);

      rt_size = CMFEM_FiniteElementSpace_GetVSize(rt_fes);
      l2_size = CMFEM_FiniteElementSpace_GetVSize(l2_fes);

      printf("Number of H(div) dofs: %d\n", rt_size);
      printf("Number of L2 dofs: %d\n", l2_size);

      CMFEM_ArrayInt_Set(&offsets, 0, 0);
      CMFEM_ArrayInt_Set(&offsets, 1, rt_size);
      CMFEM_ArrayInt_Set(&offsets, 2, rt_size + l2_size);

      x = CMFEM_Vector_ConstructSize(rt_size + l2_size);
      rhs = CMFEM_Vector_ConstructSize(rt_size + l2_size);
      delta_psi_vec = CMFEM_Vector_ConstructSize(rt_size);
      u_vec = CMFEM_Vector_ConstructSize(l2_size);
      psi_old_vec = CMFEM_Vector_ConstructSize(rt_size);
      psi_vec = CMFEM_Vector_ConstructSize(rt_size);
      u_old_vec = CMFEM_Vector_ConstructSize(l2_size);
      u_prev_newton_vec = CMFEM_Vector_ConstructSize(l2_size);
      b0_vec = CMFEM_Vector_ConstructSize(rt_size);
      b1_vec = CMFEM_Vector_ConstructSize(l2_size);
      tmp_l2_vec = CMFEM_Vector_ConstructSize(l2_size);
      tmp_l2_vec2 = CMFEM_Vector_ConstructSize(l2_size);
      tmp_u_diff = CMFEM_Vector_ConstructSize(l2_size);
      a00_diag = CMFEM_Vector_ConstructSize(rt_size);

      CMFEM_Vector_Assign(&x, 0.0);
      CMFEM_Vector_Assign(&rhs, 0.0);
      CMFEM_Vector_Assign(&delta_psi_vec, 0.0);
      CMFEM_Vector_Assign(&u_vec, 0.0);
      CMFEM_Vector_Assign(&psi_old_vec, 0.0);
      CMFEM_Vector_Assign(&psi_vec, 0.0);
      CMFEM_Vector_Assign(&u_old_vec, 0.0);

      delta_psi_gf = CMFEM_GridFunction_New(rt_fes);
      u_gf = CMFEM_GridFunction_New(l2_fes);
      psi_old_gf = CMFEM_GridFunction_New(rt_fes);
      psi_gf = CMFEM_GridFunction_New(rt_fes);
      u_old_gf = CMFEM_GridFunction_New(l2_fes);
      u_tmp_gf = CMFEM_GridFunction_New(l2_fes);

      CMFEM_GridFunction_Assign(delta_psi_gf, 0.0);
      CMFEM_GridFunction_Assign(u_gf, 0.0);
      CMFEM_GridFunction_Assign(psi_old_gf, 0.0);
      CMFEM_GridFunction_Assign(psi_gf, 0.0);
      CMFEM_GridFunction_Assign(u_old_gf, 0.0);
      CMFEM_GridFunction_Assign(u_tmp_gf, 0.0);

      // 4. Assemble the constant divergence block once.
      {
         CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
         a10 = CMFEM_MixedBilinearForm_New(rt_fes, l2_fes);
         CMFEM_MixedBilinearForm_AddDomainIntegratorVfdiCc(a10, one);
         CMFEM_MixedBilinearForm_Assemble(a10);
         CMFEM_MixedBilinearForm_Finalize(a10);
         a10_matrix = CMFEM_MixedBilinearForm_SpMat(a10);
         a01_matrix = CMFEM_TransposeSm(a10_matrix);
         CMFEM_ConstantCoefficient_Delete(one);
      }

      zero_cf = CMFEM_FunctionCoefficient_New(zero_fun, NULL);

      // 5. Run the outer proximal and inner Newton iterations.
      for (i = 0; i < max_it; i++)
      {
         int j;
         outer_iterations = i + 1;
         u_prev_newton_vec = CMFEM_Vector_Copy(&u_old_vec);

         printf("\nOUTER ITERATION %d\n", i + 1);

         for (j = 0; j < 5; j++)
         {
            CMFEM_BlockDiagonalPreconditioner *prec;
            CMFEM_DSmoother *a00_inv;
            CMFEM_GSSmoother *s_inv;
            CMFEM_BlockOperator *a_block;
            CMFEM_SparseMatrix *a00_matrix;
            CMFEM_SparseMatrix *s_matrix;
            CMFEM_LinearForm *b1_const_form;
            CMFEM_ConstantCoefficient *neg_alpha_cf;
            double newton_update_size;

            total_iterations++;

            CMFEM_AssembleIsomorphismLinearFormRtGf(rt_fes, psi_gf, &b0_vec);
            CMFEM_Vector_Scale(&b0_vec, -1.0);

            neg_alpha_cf = CMFEM_ConstantCoefficient_New(-alpha_k);
            b1_const_form = CMFEM_LinearForm_New(l2_fes);
            CMFEM_LinearForm_AddDomainIntegratorDliCc(b1_const_form, neg_alpha_cf);
            CMFEM_LinearForm_Assemble(b1_const_form);
            CMFEM_LinearForm_CopyToVector(b1_const_form, &b1_vec);
            CMFEM_LinearForm_Delete(b1_const_form);
            CMFEM_ConstantCoefficient_Delete(neg_alpha_cf);

            CMFEM_SparseMatrix_Mult(a10_matrix, &psi_old_vec, &tmp_l2_vec);
            CMFEM_SparseMatrix_Mult(a10_matrix, &psi_vec, &tmp_l2_vec2);
            CMFEM_Vector_Add(&b1_vec, &tmp_l2_vec);
            CMFEM_Vector_AddScaled(&b1_vec, -1.0, &tmp_l2_vec2);

            for (int row = 0; row < rt_size; row++)
            {
               CMFEM_Vector_Set(&rhs, row, CMFEM_Vector_Get(&b0_vec, row));
            }
            for (int row = 0; row < l2_size; row++)
            {
               CMFEM_Vector_Set(&rhs, rt_size + row, CMFEM_Vector_Get(&b1_vec, row));
            }

            a00_matrix = CMFEM_AssembleDIsomorphismMassMatrixRtGf(rt_fes, psi_gf,
                                                                  eps);
            CMFEM_SparseMatrix_GetDiag(a00_matrix, &a00_diag);
            for (int row = 0; row < rt_size; row++)
            {
               const double diag_entry = CMFEM_Vector_Get(&a00_diag, row);
               if (diag_entry != 0.0)
               {
                  CMFEM_Vector_Set(&a00_diag, row, 1.0 / diag_entry);
               }
            }
            s_matrix = CMFEM_MultAtDASmVec(a01_matrix, &a00_diag);

            a00_inv = CMFEM_DSmoother_NewSm(a00_matrix);
            s_inv = CMFEM_GSSmoother_NewSm(s_matrix);
            CMFEM_DSmoother_SetIterativeMode(a00_inv, 0);
            CMFEM_GSSmoother_SetIterativeMode(s_inv, 0);

            prec = CMFEM_BlockDiagonalPreconditioner_NewAi(&offsets);
            CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockDs(prec, 0, a00_inv);
            CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockGs(prec, 1, s_inv);

            a_block = CMFEM_BlockOperator_NewAiAi(&offsets, &offsets);
            CMFEM_BlockOperator_SetBlockSm(a_block, 0, 0, a00_matrix);
            CMFEM_BlockOperator_SetBlockSm(a_block, 1, 0, a10_matrix);
            CMFEM_BlockOperator_SetBlockSm(a_block, 0, 1, a01_matrix);

            CMFEM_MINRESBopBdp(a_block, prec, &rhs, &x, 0, 2000, 1e-12, 0.0);

            for (int row = 0; row < rt_size; row++)
            {
               CMFEM_Vector_Set(&delta_psi_vec, row, CMFEM_Vector_Get(&x, row));
            }
            for (int row = 0; row < l2_size; row++)
            {
               CMFEM_Vector_Set(&u_vec, row, CMFEM_Vector_Get(&x, rt_size + row));
            }
            CMFEM_GridFunction_SetFromVec(delta_psi_gf, &delta_psi_vec);
            CMFEM_GridFunction_SetFromVec(u_gf, &u_vec);

            tmp_u_diff = CMFEM_Vector_Copy(&u_prev_newton_vec);
            CMFEM_Vector_AddScaled(&tmp_u_diff, -1.0, &u_vec);
            CMFEM_GridFunction_SetFromVec(u_tmp_gf, &tmp_u_diff);
            newton_update_size =
               CMFEM_GridFunction_ComputeL2ErrorFc(u_tmp_gf, zero_cf);
            u_prev_newton_vec = CMFEM_Vector_Copy(&u_vec);

            CMFEM_Vector_AddScaled(&psi_vec, newton_scaling, &delta_psi_vec);
            CMFEM_GridFunction_SetFromVec(psi_gf, &psi_vec);

            if (visualization)
            {
               CMFEM_SendSolutionToGLVis(mesh, u_gf, "localhost", 19916);
            }

            printf("Newton_update_size = %.8g\n", newton_update_size);

            CMFEM_BlockOperator_Delete(a_block);
            CMFEM_BlockDiagonalPreconditioner_Delete(prec);
            CMFEM_GSSmoother_Delete(s_inv);
            CMFEM_DSmoother_Delete(a00_inv);
            CMFEM_SparseMatrix_Delete(s_matrix);
            CMFEM_SparseMatrix_Delete(a00_matrix);

            if (newton_update_size < increment_u)
            {
               break;
            }
         }

         tmp_u_diff = CMFEM_Vector_Copy(&u_vec);
         CMFEM_Vector_AddScaled(&tmp_u_diff, -1.0, &u_old_vec);
         CMFEM_GridFunction_SetFromVec(u_tmp_gf, &tmp_u_diff);
         increment_u = CMFEM_GridFunction_ComputeL2ErrorFc(u_tmp_gf, zero_cf);

         printf("Number of Newton iterations = %d\n", j + 1);
         printf("Increment (|| u_h - u_h_prvs||) = %.8g\n", increment_u);

         u_old_vec = CMFEM_Vector_Copy(&u_vec);
         psi_old_vec = CMFEM_Vector_Copy(&psi_vec);
         CMFEM_GridFunction_SetFromVec(u_old_gf, &u_old_vec);
         CMFEM_GridFunction_SetFromVec(psi_old_gf, &psi_old_vec);

         if (increment_u < tol || i == max_it - 1)
         {
            break;
         }

         alpha_k *= max_double(growth_rate, 1.0);
      }

      printf("\n Outer iterations: %d\n", outer_iterations);
      printf(" Total iterations: %d\n", total_iterations);
      printf(" Total dofs:       %d\n", rt_size + l2_size);

      CMFEM_FunctionCoefficient_Delete(zero_cf);
      CMFEM_SparseMatrix_Delete(a01_matrix);
      CMFEM_MixedBilinearForm_Delete(a10);
      CMFEM_GridFunction_Delete(u_tmp_gf);
      CMFEM_GridFunction_Delete(u_old_gf);
      CMFEM_GridFunction_Delete(psi_gf);
      CMFEM_GridFunction_Delete(psi_old_gf);
      CMFEM_GridFunction_Delete(u_gf);
      CMFEM_GridFunction_Delete(delta_psi_gf);
      CMFEM_Vector_Destroy(&a00_diag);
      CMFEM_Vector_Destroy(&tmp_u_diff);
      CMFEM_Vector_Destroy(&tmp_l2_vec2);
      CMFEM_Vector_Destroy(&tmp_l2_vec);
      CMFEM_Vector_Destroy(&b1_vec);
      CMFEM_Vector_Destroy(&b0_vec);
      CMFEM_Vector_Destroy(&u_prev_newton_vec);
      CMFEM_Vector_Destroy(&u_old_vec);
      CMFEM_Vector_Destroy(&psi_vec);
      CMFEM_Vector_Destroy(&psi_old_vec);
      CMFEM_Vector_Destroy(&u_vec);
      CMFEM_Vector_Destroy(&delta_psi_vec);
      CMFEM_Vector_Destroy(&rhs);
      CMFEM_Vector_Destroy(&x);
      CMFEM_ArrayInt_Destroy(&offsets);
      CMFEM_FiniteElementSpace_Delete(l2_fes);
      CMFEM_FiniteElementSpace_Delete(rt_fes);
      CMFEM_L2FeCollection_Delete(l2_fec);
      CMFEM_RtFeCollection_Delete(rt_fec);
      CMFEM_Mesh_Delete(mesh);
      (void)sdim;
   }

   return 0;
}
