//                              MFEM C Example 24
//
// Description: This example illustrates mixed finite element spaces and three
//              related projection problems for grad, curl, and div.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

static int dim_ = 0;
static double freq_ = 1.0;
static double kappa_ = 0.0;
static const double cmfem_pi_ = 3.14159265358979323846;

static double p_exact(const CMFEM_Vector *x, void *context)
{
   (void)context;
   if (dim_ == 3)
   {
      return sin(CMFEM_Vector_Get(x, 0)) *
             sin(CMFEM_Vector_Get(x, 1)) *
             sin(CMFEM_Vector_Get(x, 2));
   }
   if (dim_ == 2)
   {
      return sin(CMFEM_Vector_Get(x, 0)) * sin(CMFEM_Vector_Get(x, 1));
   }
   return 0.0;
}

static void gradp_exact(const CMFEM_Vector *x, CMFEM_Vector *value,
                        void *context)
{
   (void)context;
   if (dim_ == 3)
   {
      CMFEM_Vector_Set(value, 0,
                       cos(CMFEM_Vector_Get(x, 0)) *
                       sin(CMFEM_Vector_Get(x, 1)) *
                       sin(CMFEM_Vector_Get(x, 2)));
      CMFEM_Vector_Set(value, 1,
                       sin(CMFEM_Vector_Get(x, 0)) *
                       cos(CMFEM_Vector_Get(x, 1)) *
                       sin(CMFEM_Vector_Get(x, 2)));
      CMFEM_Vector_Set(value, 2,
                       sin(CMFEM_Vector_Get(x, 0)) *
                       sin(CMFEM_Vector_Get(x, 1)) *
                       cos(CMFEM_Vector_Get(x, 2)));
      return;
   }

   CMFEM_Vector_Set(value, 0,
                    cos(CMFEM_Vector_Get(x, 0)) * sin(CMFEM_Vector_Get(x, 1)));
   CMFEM_Vector_Set(value, 1,
                    sin(CMFEM_Vector_Get(x, 0)) * cos(CMFEM_Vector_Get(x, 1)));
   if (CMFEM_Vector_Size(value) == 3)
   {
      CMFEM_Vector_Set(value, 2, 0.0);
   }
}

static double div_gradp_exact(const CMFEM_Vector *x, void *context)
{
   (void)context;
   if (dim_ == 3)
   {
      return -3.0 *
             sin(CMFEM_Vector_Get(x, 0)) *
             sin(CMFEM_Vector_Get(x, 1)) *
             sin(CMFEM_Vector_Get(x, 2));
   }
   if (dim_ == 2)
   {
      return -2.0 *
             sin(CMFEM_Vector_Get(x, 0)) *
             sin(CMFEM_Vector_Get(x, 1));
   }
   return 0.0;
}

static void v_exact(const CMFEM_Vector *x, CMFEM_Vector *value, void *context)
{
   (void)context;
   if (dim_ == 3)
   {
      CMFEM_Vector_Set(value, 0, sin(kappa_ * CMFEM_Vector_Get(x, 1)));
      CMFEM_Vector_Set(value, 1, sin(kappa_ * CMFEM_Vector_Get(x, 2)));
      CMFEM_Vector_Set(value, 2, sin(kappa_ * CMFEM_Vector_Get(x, 0)));
      return;
   }

   CMFEM_Vector_Set(value, 0, sin(kappa_ * CMFEM_Vector_Get(x, 1)));
   CMFEM_Vector_Set(value, 1, sin(kappa_ * CMFEM_Vector_Get(x, 0)));
   if (CMFEM_Vector_Size(value) == 3)
   {
      CMFEM_Vector_Set(value, 2, 0.0);
   }
}

static void curlv_exact(const CMFEM_Vector *x, CMFEM_Vector *value,
                        void *context)
{
   int i;
   (void)context;
   if (dim_ == 3)
   {
      CMFEM_Vector_Set(value, 0, -kappa_ * cos(kappa_ * CMFEM_Vector_Get(x, 2)));
      CMFEM_Vector_Set(value, 1, -kappa_ * cos(kappa_ * CMFEM_Vector_Get(x, 0)));
      CMFEM_Vector_Set(value, 2, -kappa_ * cos(kappa_ * CMFEM_Vector_Get(x, 1)));
      return;
   }

   for (i = 0; i < CMFEM_Vector_Size(value); i++)
   {
      CMFEM_Vector_Set(value, i, 0.0);
   }
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("beam-hex.mesh");
   int order = 1;
   int prob = 0;
   int static_cond = 0;
   int pa = 0;
   const char *device_config = "cpu";
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

      parsed = cmfem_parse_int_option(argc, argv, &i, "-p", "--problem-type",
                                      &prob);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-sc", "--static-condensation",
                                       "-no-sc", "--no-static-condensation",
                                       &static_cond);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-pa", "--partial-assembly",
                                       "-no-pa", "--no-partial-assembly",
                                       &pa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device",
                                         &device_config);
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
   kappa_ = freq_ * cmfem_pi_;

   // 2. Enable hardware devices such as GPUs and programming models such as
   //    OpenMP based on command line options.
   CMFEM_Device *device = CMFEM_Device_New(device_config);
   CMFEM_Device_Print(device);

   // 3. Read the mesh from the given mesh file.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   dim_ = CMFEM_Mesh_Dimension(mesh);
   {
      int ref_levels = cmfem_uniform_refinement_levels(
                          50000.0,
                          CMFEM_Mesh_GetNE(mesh),
                          dim_);
      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }
   }

   // 4. Define the trial and test finite element spaces for the selected
   //    grad/curl/div projection problem.
   {
      const int sdim = CMFEM_Mesh_SpaceDimension(mesh);
      CMFEM_H1FeCollection *trial_h1 = NULL;
      CMFEM_NdFeCollection *trial_nd = NULL;
      CMFEM_RtFeCollection *trial_rt = NULL;
      CMFEM_NdFeCollection *test_nd = NULL;
      CMFEM_RtFeCollection *test_rt = NULL;
      CMFEM_L2FeCollection *test_l2 = NULL;
      CMFEM_FiniteElementSpace *trial_fes = NULL;
      CMFEM_FiniteElementSpace *test_fes = NULL;
      CMFEM_GridFunction *gftrial = NULL;
      CMFEM_GridFunction *x = NULL;
      CMFEM_GridFunction *discrete_interpolant = NULL;
      CMFEM_GridFunction *exact_proj = NULL;
      CMFEM_BilinearForm *a = NULL;
      CMFEM_MixedBilinearForm *a_mixed = NULL;
      CMFEM_DiscreteLinearOperator *dlo = NULL;
      CMFEM_ConstantCoefficient *one = NULL;
      CMFEM_FunctionCoefficient *p_coef = NULL;
      CMFEM_VectorFunctionCoefficient *gradp_coef = NULL;
      CMFEM_VectorFunctionCoefficient *v_coef = NULL;
      CMFEM_VectorFunctionCoefficient *curlv_coef = NULL;
      CMFEM_FunctionCoefficient *divgradp_coef = NULL;
      _Alignas(max_align_t) CMFEM_Vector trial_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector rhs = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector sol = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector interp = CMFEM_Vector_Construct();

      if (prob == 0)
      {
         trial_h1 = CMFEM_H1FeCollection_NewOrderDim(order, dim_);
         test_nd = CMFEM_NdFeCollection_NewOrderDim(order, dim_);
         trial_fes = CMFEM_FiniteElementSpace_NewMeshH1(mesh, trial_h1);
         test_fes = CMFEM_FiniteElementSpace_NewMeshNd(mesh, test_nd);
      }
      else if (prob == 1)
      {
         trial_nd = CMFEM_NdFeCollection_NewOrderDim(order, dim_);
         test_rt = CMFEM_RtFeCollection_NewOrderDim(order - 1, dim_);
         trial_fes = CMFEM_FiniteElementSpace_NewMeshNd(mesh, trial_nd);
         test_fes = CMFEM_FiniteElementSpace_NewMeshRt(mesh, test_rt);
      }
      else
      {
         trial_rt = CMFEM_RtFeCollection_NewOrderDim(order - 1, dim_);
         test_l2 = CMFEM_L2FeCollection_NewOrderDim(order - 1, dim_);
         trial_fes = CMFEM_FiniteElementSpace_NewMeshRt(mesh, trial_rt);
         test_fes = CMFEM_FiniteElementSpace_NewMeshL2Vdim(mesh, test_l2, 1);
      }

      if (prob == 0)
      {
         printf("Number of Nedelec finite element unknowns: %d\n",
                CMFEM_FiniteElementSpace_GetTrueVSize(test_fes));
         printf("Number of H1 finite element unknowns: %d\n",
                CMFEM_FiniteElementSpace_GetTrueVSize(trial_fes));
      }
      else if (prob == 1)
      {
         printf("Number of Nedelec finite element unknowns: %d\n",
                CMFEM_FiniteElementSpace_GetTrueVSize(trial_fes));
         printf("Number of Raviart-Thomas finite element unknowns: %d\n",
                CMFEM_FiniteElementSpace_GetTrueVSize(test_fes));
      }
      else
      {
         printf("Number of Raviart-Thomas finite element unknowns: %d\n",
                CMFEM_FiniteElementSpace_GetTrueVSize(trial_fes));
         printf("Number of L2 finite element unknowns: %d\n",
                CMFEM_FiniteElementSpace_GetTrueVSize(test_fes));
      }

      // 5. Project the exact trial-space field into the corresponding finite
      //    element grid function.
      gftrial = CMFEM_GridFunction_New(trial_fes);
      x = CMFEM_GridFunction_New(test_fes);
      discrete_interpolant = CMFEM_GridFunction_New(test_fes);
      exact_proj = CMFEM_GridFunction_New(test_fes);
      p_coef = CMFEM_FunctionCoefficient_New(p_exact, NULL);
      gradp_coef = CMFEM_VectorFunctionCoefficient_New(sdim, gradp_exact, NULL);
      v_coef = CMFEM_VectorFunctionCoefficient_New(sdim, v_exact, NULL);
      curlv_coef = CMFEM_VectorFunctionCoefficient_New(sdim, curlv_exact, NULL);
      divgradp_coef = CMFEM_FunctionCoefficient_New(div_gradp_exact, NULL);

      if (prob == 0)
      {
         CMFEM_GridFunction_ProjectCoefficientFc(gftrial, p_coef);
      }
      else if (prob == 1)
      {
         CMFEM_GridFunction_ProjectCoefficientVfc(gftrial, v_coef);
      }
      else
      {
         CMFEM_GridFunction_ProjectCoefficientVfc(gftrial, gradp_coef);
      }
      CMFEM_GridFunction_SetTrueVector(gftrial);
      CMFEM_GridFunction_SetFromTrueVector(gftrial);

      // 6. Set up the bilinear forms that define the projection solve.
      one = CMFEM_ConstantCoefficient_New(1.0);
      a = CMFEM_BilinearForm_New(test_fes);
      a_mixed = CMFEM_MixedBilinearForm_New(trial_fes, test_fes);
      if (pa)
      {
         CMFEM_BilinearForm_SetAssemblyLevelPartial(a);
         CMFEM_MixedBilinearForm_SetAssemblyLevelPartial(a_mixed);
      }

      if (prob == 0)
      {
         CMFEM_BilinearForm_AddDomainIntegratorVmi(a, one);
         CMFEM_MixedBilinearForm_AddDomainIntegratorMvgiCc(a_mixed, one);
      }
      else if (prob == 1)
      {
         CMFEM_BilinearForm_AddDomainIntegratorVmi(a, one);
         CMFEM_MixedBilinearForm_AddDomainIntegratorMvciCc(a_mixed, one);
      }
      else
      {
         CMFEM_BilinearForm_AddDomainIntegratorMiCc(a, one);
         CMFEM_MixedBilinearForm_AddDomainIntegratorVfdiCc(a_mixed, one);
      }
      if (static_cond)
      {
         CMFEM_BilinearForm_EnableStaticCondensation(a);
      }

      CMFEM_BilinearForm_Assemble(a);
      if (!pa)
      {
         CMFEM_BilinearForm_Finalize(a);
      }
      CMFEM_MixedBilinearForm_Assemble(a_mixed);
      if (!pa)
      {
         CMFEM_MixedBilinearForm_Finalize(a_mixed);
      }

      // 7. Assemble the mixed right-hand side and solve the projection system
      //    with CG and Jacobi preconditioning.
      trial_vec = CMFEM_Vector_ConstructSize(CMFEM_FiniteElementSpace_GetVSize(
                                                trial_fes));
      rhs = CMFEM_Vector_ConstructSize(CMFEM_FiniteElementSpace_GetVSize(test_fes));
      sol = CMFEM_Vector_ConstructSize(CMFEM_FiniteElementSpace_GetVSize(test_fes));
      interp = CMFEM_Vector_ConstructSize(CMFEM_FiniteElementSpace_GetVSize(
                                             test_fes));
      CMFEM_GridFunction_CopyToVec(gftrial, &trial_vec);
      CMFEM_Vector_Assign(&rhs, 0.0);
      CMFEM_Vector_Assign(&sol, 0.0);
      CMFEM_Vector_Assign(&interp, 0.0);
      CMFEM_MixedBilinearForm_Mult(a_mixed, &trial_vec, &rhs);

      if (pa)
      {
         _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list =
            CMFEM_ArrayInt_Construct();
         CMFEM_OperatorJacobiSmoother *jacobi =
            CMFEM_OperatorJacobiSmoother_NewBf(a, &ess_tdof_list);
         CMFEM_PCGBfOjs(a, jacobi, &rhs, &sol, 1, 1000, 1e-12, 0.0);
         CMFEM_OperatorJacobiSmoother_Delete(jacobi);
         CMFEM_ArrayInt_Destroy(&ess_tdof_list);
      }
      else
      {
         CMFEM_DSmoother *jacobi = CMFEM_DSmoother_NewBf(a);
         CMFEM_PCGBfDs(a, jacobi, &rhs, &sol, 1, 1000, 1e-12, 0.0);
         CMFEM_DSmoother_Delete(jacobi);
      }
      CMFEM_GridFunction_SetFromVec(x, &sol);

      // 8. Apply the discrete interpolator that represents grad, curl, or div
      //    directly between the finite element spaces.
      dlo = CMFEM_DiscreteLinearOperator_New(trial_fes, test_fes);
      if (prob == 0)
      {
         CMFEM_DiscreteLinearOperator_AddDomainInterpolatorGi(dlo);
      }
      else if (prob == 1)
      {
         CMFEM_DiscreteLinearOperator_AddDomainInterpolatorCi(dlo);
      }
      else
      {
         CMFEM_DiscreteLinearOperator_AddDomainInterpolatorDi(dlo);
      }
      CMFEM_DiscreteLinearOperator_Assemble(dlo);
      CMFEM_DiscreteLinearOperator_Mult(dlo, &trial_vec, &interp);
      CMFEM_GridFunction_SetFromVec(discrete_interpolant, &interp);

      // 9. Project the exact target field into the test space.
      if (prob == 0)
      {
         CMFEM_GridFunction_ProjectCoefficientVfc(exact_proj, gradp_coef);
      }
      else if (prob == 1)
      {
         CMFEM_GridFunction_ProjectCoefficientVfc(exact_proj, curlv_coef);
      }
      else
      {
         CMFEM_GridFunction_ProjectCoefficientFc(exact_proj, divgradp_coef);
      }
      CMFEM_GridFunction_SetTrueVector(exact_proj);
      CMFEM_GridFunction_SetFromTrueVector(exact_proj);

      // 10. Compute and print the L2 error for the solve, the interpolant, and
      //     the projection of the exact target field.
      if (prob == 0)
      {
         const double err_sol = CMFEM_GridFunction_ComputeL2ErrorVfc(x, gradp_coef);
         const double err_interp = CMFEM_GridFunction_ComputeL2ErrorVfc(
                                      discrete_interpolant,
                                      gradp_coef);
         const double err_proj = CMFEM_GridFunction_ComputeL2ErrorVfc(
                                    exact_proj,
                                    gradp_coef);
         printf("\n Solution of (E_h,v) = (grad p_h,v) for E_h and v in H(curl): || E_h - grad p ||_{L_2} = %.14g\n\n",
                err_sol);
         printf(" Gradient interpolant E_h = grad p_h in H(curl): || E_h - grad p ||_{L_2} = %.14g\n\n",
                err_interp);
         printf(" Projection E_h of exact grad p in H(curl): || E_h - grad p ||_{L_2} = %.14g\n\n",
                err_proj);
      }
      else if (prob == 1)
      {
         const double err_sol = CMFEM_GridFunction_ComputeL2ErrorVfc(x, curlv_coef);
         const double err_interp = CMFEM_GridFunction_ComputeL2ErrorVfc(
                                      discrete_interpolant,
                                      curlv_coef);
         const double err_proj = CMFEM_GridFunction_ComputeL2ErrorVfc(
                                    exact_proj,
                                    curlv_coef);
         printf("\n Solution of (E_h,w) = (curl v_h,w) for E_h and w in H(div): || E_h - curl v ||_{L_2} = %.14g\n\n",
                err_sol);
         printf(" Curl interpolant E_h = curl v_h in H(div): || E_h - curl v ||_{L_2} = %.14g\n\n",
                err_interp);
         printf(" Projection E_h of exact curl v in H(div): || E_h - curl v ||_{L_2} = %.14g\n\n",
                err_proj);
      }
      else
      {
         const int quadrature_order = order > 0 ? 2 * order + 1 : 2;
         const double err_sol = CMFEM_GridFunction_ComputeL2ErrorFcOrder(
                                   x,
                                   divgradp_coef,
                                   quadrature_order);
         const double err_interp = CMFEM_GridFunction_ComputeL2ErrorFcOrder(
                                      discrete_interpolant,
                                      divgradp_coef,
                                      quadrature_order);
         const double err_proj = CMFEM_GridFunction_ComputeL2ErrorFcOrder(
                                    exact_proj,
                                    divgradp_coef,
                                    quadrature_order);
         printf("\n Solution of (f_h,q) = (div v_h,q) for f_h and q in L_2: || f_h - div v ||_{L_2} = %.14g\n\n",
                err_sol);
         printf(" Divergence interpolant f_h = div v_h in L_2: || f_h - div v ||_{L_2} = %.14g\n\n",
                err_interp);
         printf(" Projection f_h of exact div v in L_2: || f_h - div v ||_{L_2} = %.14g\n\n",
                err_proj);
      }

      // 11. Save the refined mesh and the projected solution.
      CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
      CMFEM_GridFunction_Save(x, "sol.gf", 8);

      // 12. Send the projected solution to a GLVis server.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
      }

      CMFEM_Vector_Destroy(&interp);
      CMFEM_Vector_Destroy(&sol);
      CMFEM_Vector_Destroy(&rhs);
      CMFEM_Vector_Destroy(&trial_vec);
      CMFEM_DiscreteLinearOperator_Delete(dlo);
      CMFEM_MixedBilinearForm_Delete(a_mixed);
      CMFEM_BilinearForm_Delete(a);
      CMFEM_GridFunction_Delete(exact_proj);
      CMFEM_GridFunction_Delete(discrete_interpolant);
      CMFEM_GridFunction_Delete(x);
      CMFEM_GridFunction_Delete(gftrial);
      CMFEM_FunctionCoefficient_Delete(divgradp_coef);
      CMFEM_VectorFunctionCoefficient_Delete(curlv_coef);
      CMFEM_VectorFunctionCoefficient_Delete(v_coef);
      CMFEM_VectorFunctionCoefficient_Delete(gradp_coef);
      CMFEM_FunctionCoefficient_Delete(p_coef);
      CMFEM_ConstantCoefficient_Delete(one);
      CMFEM_FiniteElementSpace_Delete(test_fes);
      CMFEM_FiniteElementSpace_Delete(trial_fes);
      CMFEM_L2FeCollection_Delete(test_l2);
      CMFEM_RtFeCollection_Delete(test_rt);
      CMFEM_NdFeCollection_Delete(test_nd);
      CMFEM_RtFeCollection_Delete(trial_rt);
      CMFEM_NdFeCollection_Delete(trial_nd);
      CMFEM_H1FeCollection_Delete(trial_h1);
   }

   CMFEM_Mesh_Delete(mesh);
   CMFEM_Device_Delete(device);
   return 0;
}
