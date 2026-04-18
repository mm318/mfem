//                               MFEM C Example 16
//
// Description: This example solves a nonlinear heat equation in time by
//              building a C-side conduction operator and advancing it with an
//              MFEM ODE solver through the CMFEM wrapper layer.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

typedef struct ConductionOperatorContext
{
   CMFEM_FiniteElementSpace *fespace;
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list;
   _Alignas(max_align_t) CMFEM_SparseMatrix mmat;
   _Alignas(max_align_t) CMFEM_SparseMatrix kmat;
   _Alignas(max_align_t) CMFEM_Vector z;
   CMFEM_BilinearForm *m_form;
   CMFEM_BilinearForm *k_form;
   CMFEM_SparseMatrix *tmat;
   CMFEM_DSmoother *m_prec;
   CMFEM_DSmoother *t_prec;
   double alpha;
   double kappa;
   double current_gamma;
   int solve_implicit_state;
} ConductionOperatorContext;

static double initial_temperature(const CMFEM_Vector *x, void *context)
{
   int i;
   double norm2 = 0.0;
   (void)context;

   for (i = 0; i < CMFEM_Vector_Size(x); i++)
   {
      const double value = CMFEM_Vector_Get(x, i);
      norm2 += value * value;
   }

   return sqrt(norm2) < 0.5 ? 2.0 : 1.0;
}

static void conduction_reset_implicit(ConductionOperatorContext *oper)
{
   if (oper->t_prec != NULL)
   {
      CMFEM_DSmoother_Delete(oper->t_prec);
      oper->t_prec = NULL;
   }
   if (oper->tmat != NULL)
   {
      CMFEM_SparseMatrix_Delete(oper->tmat);
      oper->tmat = NULL;
   }
}

static void conduction_set_parameters(ConductionOperatorContext *oper,
                                      const CMFEM_Vector *u)
{
   CMFEM_GridFunction *u_alpha_gf = CMFEM_GridFunction_New(oper->fespace);
   CMFEM_GridFunctionCoefficient *u_coeff = NULL;
   int i;

   CMFEM_GridFunction_SetFromTrueDofs(u_alpha_gf, u);
   for (i = 0; i < CMFEM_GridFunction_Size(u_alpha_gf); i++)
   {
      const double value = oper->kappa +
                           oper->alpha * CMFEM_GridFunction_Get(u_alpha_gf, i);
      CMFEM_GridFunction_Set(u_alpha_gf, i, value);
   }

   if (oper->k_form != NULL)
   {
      CMFEM_BilinearForm_Delete(oper->k_form);
   }
   oper->k_form = CMFEM_BilinearForm_New(oper->fespace);
   u_coeff = CMFEM_GridFunctionCoefficient_NewGf(u_alpha_gf);
   CMFEM_BilinearForm_AddDomainIntegratorDiGfc(oper->k_form, u_coeff);
   CMFEM_BilinearForm_Assemble(oper->k_form);
   CMFEM_BilinearForm_FormSystemMatrixSm(oper->k_form, &oper->ess_tdof_list,
                                         &oper->kmat);

   CMFEM_GridFunctionCoefficient_Delete(u_coeff);
   CMFEM_GridFunction_Delete(u_alpha_gf);

   conduction_reset_implicit(oper);
}

static void conduction_init(ConductionOperatorContext *oper,
                            CMFEM_FiniteElementSpace *fespace,
                            double alpha,
                            double kappa,
                            int solve_implicit_state,
                            const CMFEM_Vector *u)
{
   CMFEM_ConstantCoefficient *one = NULL;

   oper->fespace = fespace;
   oper->ess_tdof_list = CMFEM_ArrayInt_Construct();
   oper->mmat = CMFEM_SparseMatrix_Construct();
   oper->kmat = CMFEM_SparseMatrix_Construct();
   oper->z = CMFEM_Vector_ConstructSize(
                CMFEM_FiniteElementSpace_GetTrueVSize(fespace));
   oper->m_form = NULL;
   oper->k_form = NULL;
   oper->tmat = NULL;
   oper->m_prec = NULL;
   oper->t_prec = NULL;
   oper->alpha = alpha;
   oper->kappa = kappa;
   oper->current_gamma = -1.0;
   oper->solve_implicit_state = solve_implicit_state;

   oper->m_form = CMFEM_BilinearForm_New(fespace);
   one = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_BilinearForm_AddDomainIntegratorMiCc(oper->m_form, one);
   CMFEM_BilinearForm_Assemble(oper->m_form);
   CMFEM_BilinearForm_FormSystemMatrixSm(oper->m_form, &oper->ess_tdof_list,
                                         &oper->mmat);
   CMFEM_ConstantCoefficient_Delete(one);

   oper->m_prec = CMFEM_DSmoother_NewSm(&oper->mmat);

   conduction_set_parameters(oper, u);
}

static void conduction_destroy(ConductionOperatorContext *oper)
{
   conduction_reset_implicit(oper);
   if (oper->m_prec != NULL)
   {
      CMFEM_DSmoother_Delete(oper->m_prec);
   }
   if (oper->k_form != NULL)
   {
      CMFEM_BilinearForm_Delete(oper->k_form);
   }
   if (oper->m_form != NULL)
   {
      CMFEM_BilinearForm_Delete(oper->m_form);
   }
   CMFEM_Vector_Destroy(&oper->z);
   CMFEM_SparseMatrix_Destroy(&oper->kmat);
   CMFEM_SparseMatrix_Destroy(&oper->mmat);
   CMFEM_ArrayInt_Destroy(&oper->ess_tdof_list);
}

static void conduction_mult(const CMFEM_Vector *u,
                            CMFEM_Vector *du_dt,
                            double time,
                            void *context)
{
   ConductionOperatorContext *oper = (ConductionOperatorContext *)context;
   (void)time;

   // Compute du_dt = M^{-1}(-K u).
   CMFEM_SparseMatrix_Mult(&oper->kmat, u, &oper->z);
   CMFEM_Vector_Neg(&oper->z);
   CMFEM_PCGSmDs(&oper->mmat,
                 oper->m_prec,
                 &oper->z,
                 du_dt,
                 0,
                 30,
                 1.0e-8,
                 0.0);
}

static void conduction_implicit_solve(double gamma,
                                      const CMFEM_Vector *u,
                                      CMFEM_Vector *k,
                                      double time,
                                      void *context)
{
   ConductionOperatorContext *oper = (ConductionOperatorContext *)context;
   const double tolerance = 1.0e-14 * (1.0 + fabs(gamma));
   (void)time;

   // Solve either for the stage slope k = du/dt or for the stage state k = u_s.
   if (oper->tmat == NULL || fabs(gamma - oper->current_gamma) > tolerance)
   {
      conduction_reset_implicit(oper);
      oper->tmat = CMFEM_AddSmSm(1.0, &oper->mmat, gamma, &oper->kmat);
      oper->t_prec = CMFEM_DSmoother_NewSm(oper->tmat);
      oper->current_gamma = gamma;
   }

   if (oper->solve_implicit_state)
   {
      CMFEM_SparseMatrix_Mult(&oper->mmat, u, &oper->z);
   }
   else
   {
      CMFEM_SparseMatrix_Mult(&oper->kmat, u, &oper->z);
      CMFEM_Vector_Neg(&oper->z);
   }

   CMFEM_PCGSmDs(oper->tmat,
                 oper->t_prec,
                 &oper->z,
                 k,
                 0,
                 100,
                 1.0e-8,
                 0.0);
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int ref_levels = 2;
   int order = 2;
   int ode_solver_type = 23;
   double t_final = 0.5;
   double dt = 1.0e-2;
   double alpha = 1.0e-2;
   double kappa = 0.5;
   int visualization = 1;
   int visit = 0;
   int vis_steps = 5;
   int solve_implicit_state = 0;
   const int precision = 8;
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

      parsed = cmfem_parse_int_option(argc, argv, &i, "-s", "--ode-solver",
                                      &ode_solver_type);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-tf", "--t-final",
                                         &t_final);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-dt", "--time-step",
                                         &dt);
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
                                       "-imp-state", "--implicit-state",
                                       "-imp-slope", "--implicit-slope",
                                       &solve_implicit_state);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-visit", "--visit-datafiles",
                                       "-no-visit", "--no-visit-datafiles",
                                       &visit);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-vs",
                                      "--visualization-steps", &vis_steps);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   if (visit)
   {
      fprintf(stderr,
              "VisIt output is currently disabled in examples/c/ex16.c to avoid instability in the C port.\n");
      visit = 0;
   }

   // 2. Read and refine the mesh.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   const int dim = CMFEM_Mesh_Dimension(mesh);
   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }

   // 3. Define the ODE solver used for time integration.
   CMFEM_ODESolver *ode_solver = CMFEM_ODESolver_NewType(ode_solver_type);
   if (ode_solver == NULL)
   {
      fprintf(stderr, "Unsupported ODE solver type: %d\n", ode_solver_type);
      CMFEM_Mesh_Delete(mesh);
      return 1;
   }

   // 4. Define the finite element space and the initial temperature.
   CMFEM_H1FeCollection *fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh,
                                                                          fec);
   CMFEM_GridFunction *u_gf = CMFEM_GridFunction_New(fespace);
   CMFEM_FunctionCoefficient *u0 =
      CMFEM_FunctionCoefficient_New(initial_temperature, NULL);
   _Alignas(max_align_t) CMFEM_Vector u = CMFEM_Vector_Construct();
   ConductionOperatorContext oper;
   CMFEM_TimeDependentOperator *oper_adapter = NULL;
   CMFEM_VisItDataCollection *visit_dc = NULL;

   printf("Number of temperature unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   CMFEM_GridFunction_ProjectCoefficientFc(u_gf, u0);
   CMFEM_GridFunction_GetTrueDofs(u_gf, &u);

   // 5. Initialize the conduction operator and choose the implicit stage mode.
   conduction_init(&oper, fespace, alpha, kappa, solve_implicit_state, &u);
   oper_adapter = CMFEM_TimeDependentOperator_New(
                     CMFEM_FiniteElementSpace_GetTrueVSize(fespace),
                     CMFEM_TIME_DEPENDENT_OPERATOR_IMPLICIT,
                     conduction_mult,
                     &oper,
                     conduction_implicit_solve,
                     &oper);
   CMFEM_TimeDependentOperator_SetImplicitVariableType(
      oper_adapter,
      solve_implicit_state ? CMFEM_IMPLICIT_VARIABLE_STATE
      : CMFEM_IMPLICIT_VARIABLE_SLOPE);

   // 6. Save the initial solution and initialize visualization output.
   CMFEM_GridFunction_SetFromTrueDofs(u_gf, &u);
   CMFEM_Mesh_Print(mesh, "ex16.mesh", precision);
   CMFEM_GridFunction_Save(u_gf, "ex16-init.gf", precision);

   if (visit)
   {
      visit_dc = CMFEM_VisItDataCollection_New("Example16", mesh);
      CMFEM_VisItDataCollection_RegisterFieldGf(visit_dc, "temperature", u_gf);
      CMFEM_VisItDataCollection_SetCycle(visit_dc, 0);
      CMFEM_VisItDataCollection_SetTime(visit_dc, 0.0);
      CMFEM_VisItDataCollection_Save(visit_dc);
   }

   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, u_gf, "localhost", 19916);
   }

   // 7. Perform time integration.
   {
      double t = 0.0;
      int ti = 0;
      int last_step = 0;

      CMFEM_ODESolver_Init(ode_solver, oper_adapter);

      while (!last_step)
      {
         ti++;
         if (t + dt >= t_final - dt / 2.0)
         {
            last_step = 1;
         }

         CMFEM_ODESolver_Step(ode_solver, &u, &t, &dt);

         if (last_step || (ti % vis_steps) == 0)
         {
            printf("step %d, t = %.17g\n", ti, t);

            CMFEM_GridFunction_SetFromTrueDofs(u_gf, &u);
            if (visualization)
            {
               CMFEM_SendSolutionToGLVis(mesh, u_gf, "localhost", 19916);
            }

            if (visit)
            {
               CMFEM_VisItDataCollection_SetCycle(visit_dc, ti);
               CMFEM_VisItDataCollection_SetTime(visit_dc, t);
               CMFEM_VisItDataCollection_Save(visit_dc);
            }
         }

         conduction_set_parameters(&oper, &u);
      }
   }

   // 8. Save the final solution for later visualization.
   CMFEM_GridFunction_SetFromTrueDofs(u_gf, &u);
   CMFEM_GridFunction_Save(u_gf, "ex16-final.gf", precision);

   // 9. Free the used memory.
   if (visit_dc != NULL)
   {
      CMFEM_VisItDataCollection_Delete(visit_dc);
   }
   CMFEM_TimeDependentOperator_Delete(oper_adapter);
   conduction_destroy(&oper);
   CMFEM_Vector_Destroy(&u);
   CMFEM_FunctionCoefficient_Delete(u0);
   CMFEM_GridFunction_Delete(u_gf);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_H1FeCollection_Delete(fec);
   CMFEM_ODESolver_Delete(ode_solver);
   CMFEM_Mesh_Delete(mesh);

   return 0;
}
