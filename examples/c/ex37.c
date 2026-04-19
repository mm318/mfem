//                              MFEM C Example 37
//
// Description: This example demonstrates density-filtered topology
//              optimization with a reduced-space projected mirror-descent
//              update.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

static double inv_sigmoid(double x)
{
   const double tol = 1.0e-12;
   if (x < tol) { x = tol; }
   if (x > 1.0 - tol) { x = 1.0 - tol; }
   return log(x / (1.0 - x));
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int ref_levels = 5;
   int order = 2;
   double alpha = 1.0;
   double epsilon = 0.01;
   double vol_fraction = 0.5;
   int max_it = 1000;
   double itol = 1.0e-1;
   double ntol = 1.0e-4;
   double rho_min = 1.0e-6;
   double lambda = 1.0;
   double mu = 1.0;
   int visualization = 1;
   int paraview_output = 0;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_int_option(argc, argv, &i, "-r", "--refine",
                                          &ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i,
                                         "-alpha", "--alpha-step-length",
                                         &alpha);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i,
                                         "-epsilon", "--epsilon-thickness",
                                         &epsilon);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-mi", "--max-it",
                                      &max_it);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-ntol", "--rel-tol",
                                         &ntol);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-itol", "--abs-tol",
                                         &itol);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i,
                                         "-vf", "--volume-fraction",
                                         &vol_fraction);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i,
                                         "-lambda", "--lambda",
                                         &lambda);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-mu", "--mu", &mu);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i,
                                         "-rmin", "--psi-min",
                                         &rho_min);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-pv", "--paraview",
                                       "-no-pv", "--no-paraview",
                                       &paraview_output);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Construct the beam mesh and set up the finite element spaces.
   {
      CMFEM_Mesh *mesh = CMFEM_Ex37_NewBeamMesh();
      const int dim = CMFEM_Mesh_Dimension(mesh);
      CMFEM_H1FeCollection *state_fec = NULL;
      CMFEM_H1FeCollection *filter_fec = NULL;
      CMFEM_L2FeCollection *control_fec = NULL;
      CMFEM_FiniteElementSpace *state_fes = NULL;
      CMFEM_FiniteElementSpace *filter_fes = NULL;
      CMFEM_FiniteElementSpace *control_fes = NULL;
      CMFEM_GridFunction *u = NULL;
      CMFEM_GridFunction *psi = NULL;
      CMFEM_GridFunction *psi_old = NULL;
      CMFEM_GridFunction *rho_filter = NULL;
      CMFEM_GridFunction *rho_gf = NULL;
      CMFEM_GridFunction *grad = NULL;
      CMFEM_GridFunction *w_filter = NULL;
      CMFEM_GridFunction *onegf = NULL;
      CMFEM_GridFunction *r_gf = NULL;
      CMFEM_LinearForm *vol_form = NULL;
      CMFEM_BilinearForm *mass = NULL;
      CMFEM_ConstantCoefficient *one = NULL;
      CMFEM_ParaViewDataCollection *paraview_dc = NULL;
      void *filter_solver = NULL;
      void *elasticity_solver = NULL;
      double center[2] = { 2.9, 0.5 };
      double force[2] = { 0.0, -1.0 };
      double domain_volume;
      double target_volume;
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr_filter =
         CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt empty = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_SparseMatrix M = CMFEM_SparseMatrix_Construct();
      _Alignas(max_align_t) CMFEM_Vector w_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector grad_vec = CMFEM_Vector_Construct();

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      state_fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      filter_fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      control_fec = CMFEM_L2FeCollection_NewOrderDimGl(order - 1, dim);
      state_fes = CMFEM_FiniteElementSpace_NewMeshH1VDim(mesh, state_fec, dim);
      filter_fes = CMFEM_FiniteElementSpace_NewMeshH1(mesh, filter_fec);
      control_fes = CMFEM_FiniteElementSpace_NewMeshL2Vdim(mesh, control_fec, 1);

      printf("Number of state unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(state_fes));
      printf("Number of filter unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(filter_fes));
      printf("Number of control unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(control_fes));

      u = CMFEM_GridFunction_New(state_fes);
      psi = CMFEM_GridFunction_New(control_fes);
      psi_old = CMFEM_GridFunction_New(control_fes);
      rho_filter = CMFEM_GridFunction_New(filter_fes);
      rho_gf = CMFEM_GridFunction_New(control_fes);
      grad = CMFEM_GridFunction_New(control_fes);
      w_filter = CMFEM_GridFunction_New(filter_fes);
      onegf = CMFEM_GridFunction_New(control_fes);
      r_gf = CMFEM_GridFunction_New(filter_fes);
      CMFEM_GridFunction_Assign(u, 0.0);
      CMFEM_GridFunction_Assign(rho_filter, vol_fraction);
      CMFEM_GridFunction_Assign(psi, inv_sigmoid(vol_fraction));
      CMFEM_GridFunction_Assign(psi_old, inv_sigmoid(vol_fraction));
      CMFEM_GridFunction_Assign(grad, 0.0);
      CMFEM_GridFunction_Assign(w_filter, 0.0);
      CMFEM_GridFunction_Assign(onegf, 1.0);
      CMFEM_GridFunction_Assign(r_gf, 0.0);

      if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
      {
         CMFEM_ArrayInt_Destroy(&ess_bdr);
         ess_bdr = CMFEM_ArrayInt_ConstructSize(
                      CMFEM_Mesh_BoundaryAttributesMax(mesh));
         CMFEM_ArrayInt_Assign(&ess_bdr, 0);
         CMFEM_ArrayInt_Set(&ess_bdr, 0, 1);

         CMFEM_ArrayInt_Destroy(&ess_bdr_filter);
         ess_bdr_filter = CMFEM_ArrayInt_ConstructSize(
                             CMFEM_Mesh_BoundaryAttributesMax(mesh));
         CMFEM_ArrayInt_Assign(&ess_bdr_filter, 0);
      }

      filter_solver = CMFEM_Ex37FilterSolver_NewMeshOrderEpsAi(
                         mesh,
                         order,
                         epsilon * epsilon,
                         &ess_bdr_filter);
      elasticity_solver = CMFEM_Ex37ElasticitySolver_NewMeshOrderAiCircle(
                             mesh,
                             order,
                             &ess_bdr,
                             0.05,
                             center,
                             force);

      one = CMFEM_ConstantCoefficient_New(1.0);
      mass = CMFEM_BilinearForm_New(control_fes);
      CMFEM_BilinearForm_AddDomainIntegratorIiDiMiCc(mass, one);
      CMFEM_BilinearForm_Assemble(mass);
      CMFEM_BilinearForm_FormSystemMatrixSm(mass, &empty, &M);
      CMFEM_Vector_Destroy(&w_vec);
      w_vec = CMFEM_Vector_ConstructSize(CMFEM_SparseMatrix_Height(&M));
      CMFEM_Vector_Destroy(&grad_vec);
      grad_vec = CMFEM_Vector_ConstructSize(CMFEM_SparseMatrix_Height(&M));

      vol_form = CMFEM_LinearForm_New(control_fes);
      CMFEM_LinearForm_AddDomainIntegratorDliCc(vol_form, one);
      CMFEM_LinearForm_Assemble(vol_form);
      domain_volume = CMFEM_LinearForm_EvalGf(vol_form, onegf);
      target_volume = domain_volume * vol_fraction;

      if (paraview_output)
      {
         paraview_dc = CMFEM_ParaViewDataCollection_New("ex37", mesh);
         CMFEM_Ex37_ProjectSigmoidGf(psi, rho_gf);
         CMFEM_ParaViewDataCollection_SetPrefixPath(paraview_dc, "ParaView");
         CMFEM_ParaViewDataCollection_SetLevelsOfDetail(paraview_dc, order);
         CMFEM_ParaViewDataCollection_SetDataFormatBinary(paraview_dc);
         CMFEM_ParaViewDataCollection_SetHighOrderOutput(paraview_dc, 1);
         CMFEM_ParaViewDataCollection_SetCycle(paraview_dc, 0);
         CMFEM_ParaViewDataCollection_SetTime(paraview_dc, 0.0);
         CMFEM_ParaViewDataCollection_RegisterFieldGf(paraview_dc,
                                                      "displacement",
                                                      u);
         CMFEM_ParaViewDataCollection_RegisterFieldGf(paraview_dc,
                                                      "density",
                                                      rho_gf);
         CMFEM_ParaViewDataCollection_RegisterFieldGf(paraview_dc,
                                                      "filtered_density",
                                                      rho_filter);
         CMFEM_ParaViewDataCollection_Save(paraview_dc);
      }

      // 3. Run the projected mirror-descent optimization loop.
      for (i = 1; i <= max_it; i++)
      {
         double material_volume;
         double norm_increment;
         double norm_reduced_gradient;
         double compliance;
         CMFEM_GridFunctionCoefficient *w_cf = NULL;
         CMFEM_LinearForm *w_rhs = NULL;

         if (i > 1)
         {
            alpha *= ((double)i) / ((double)i - 1.0);
         }

         printf("\nStep = %d\n", i);

         // Step 1. Solve the Helmholtz filter equation.
         CMFEM_Ex37FilterSolver_SolveSigmoidGf(filter_solver, psi, rho_filter);

         // Step 2. Solve the elasticity state equation with SIMP material law.
         CMFEM_Ex37ElasticitySolver_SolveGf(elasticity_solver,
                                            lambda,
                                            mu,
                                            rho_filter,
                                            rho_min,
                                            u);

         // Step 3. Solve the adjoint-filter equation for the reduced gradient.
         CMFEM_Ex37FilterSolver_SolveSedGfGf(filter_solver,
                                             lambda,
                                             mu,
                                             u,
                                             rho_filter,
                                             rho_min,
                                             w_filter);

         // Step 4. Project the filtered gradient back to the control space.
         w_cf = CMFEM_GridFunctionCoefficient_NewGf(w_filter);
         w_rhs = CMFEM_LinearForm_New(control_fes);
         CMFEM_LinearForm_AddDomainIntegratorDliGfc(w_rhs, w_cf);
         CMFEM_LinearForm_Assemble(w_rhs);
         CMFEM_LinearForm_CopyToVector(w_rhs, &w_vec);
         CMFEM_SparseMatrix_Mult(&M, &w_vec, &grad_vec);
         CMFEM_GridFunction_SetFromVec(grad, &grad_vec);
         CMFEM_LinearForm_Delete(w_rhs);
         CMFEM_GridFunctionCoefficient_Delete(w_cf);

         // Step 5. Apply the projected mirror-descent update.
         CMFEM_GridFunction_AddScaledGf(psi, -alpha, grad);
         material_volume = CMFEM_Ex37_ProjectSigmoidVolumeGf(psi,
                                                             target_volume,
                                                             1.0e-12,
                                                             10);

         norm_increment = CMFEM_Ex37_SigmoidDiffL1GfGf(psi, psi_old);
         norm_reduced_gradient = norm_increment / alpha;
         CMFEM_GridFunction_CopyFromGf(psi_old, psi);
         compliance = CMFEM_Ex37ElasticitySolver_ComplianceGf(elasticity_solver,
                                                              u);

         printf("norm of the reduced gradient = %.8g\n", norm_reduced_gradient);
         printf("norm of the increment = %.8g\n", norm_increment);
         printf("compliance = %.8g\n", compliance);
         printf("volume fraction = %.8g\n", material_volume / domain_volume);

         if (visualization)
         {
            CMFEM_Ex37_ProjectSimpGf(rho_filter, rho_min, r_gf);
            CMFEM_SendSolutionToGLVis(mesh, r_gf, "localhost", 19916);
         }

         if (paraview_output)
         {
            CMFEM_Ex37_ProjectSigmoidGf(psi, rho_gf);
            CMFEM_ParaViewDataCollection_SetCycle(paraview_dc, i);
            CMFEM_ParaViewDataCollection_SetTime(paraview_dc, (double)i);
            CMFEM_ParaViewDataCollection_Save(paraview_dc);
         }

         if (norm_reduced_gradient < ntol && norm_increment < itol)
         {
            break;
         }
      }

      if (paraview_dc != NULL)
      {
         CMFEM_ParaViewDataCollection_Delete(paraview_dc);
      }
      CMFEM_Ex37ElasticitySolver_Delete(elasticity_solver);
      CMFEM_Ex37FilterSolver_Delete(filter_solver);
      CMFEM_LinearForm_Delete(vol_form);
      CMFEM_BilinearForm_Delete(mass);
      CMFEM_ConstantCoefficient_Delete(one);
      CMFEM_GridFunction_Delete(r_gf);
      CMFEM_GridFunction_Delete(onegf);
      CMFEM_GridFunction_Delete(w_filter);
      CMFEM_GridFunction_Delete(grad);
      CMFEM_GridFunction_Delete(rho_gf);
      CMFEM_GridFunction_Delete(rho_filter);
      CMFEM_GridFunction_Delete(psi_old);
      CMFEM_GridFunction_Delete(psi);
      CMFEM_GridFunction_Delete(u);
      CMFEM_FiniteElementSpace_Delete(control_fes);
      CMFEM_FiniteElementSpace_Delete(filter_fes);
      CMFEM_FiniteElementSpace_Delete(state_fes);
      CMFEM_L2FeCollection_Delete(control_fec);
      CMFEM_H1FeCollection_Delete(filter_fec);
      CMFEM_H1FeCollection_Delete(state_fec);
      CMFEM_Mesh_Delete(mesh);
      CMFEM_SparseMatrix_Destroy(&M);
      CMFEM_Vector_Destroy(&grad_vec);
      CMFEM_Vector_Destroy(&w_vec);
      CMFEM_ArrayInt_Destroy(&empty);
      CMFEM_ArrayInt_Destroy(&ess_bdr_filter);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
   }

   return 0;
}
