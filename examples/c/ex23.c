//                               MFEM C Example 23
//
// Description: This example solves the wave equation with second-order time
//              integration. It builds a C wave-operator context on top of the
//              CMFEM wrappers and advances the solution with a selectable
//              second-order ODE solver.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct WaveOperatorContext
{
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list;
   _Alignas(max_align_t) CMFEM_SparseMatrix mmat;
   _Alignas(max_align_t) CMFEM_SparseMatrix kmat;
   _Alignas(max_align_t) CMFEM_Vector z;
   CMFEM_BilinearForm *m_form;
   CMFEM_BilinearForm *k_form;
   CMFEM_SparseMatrix *tmat;
   CMFEM_DSmoother *m_prec;
   CMFEM_DSmoother *t_prec;
   CMFEM_ConstantCoefficient *c2;
} WaveOperatorContext;

static double initial_solution(const CMFEM_Vector *x, void *context)
{
   int i;
   double norm2 = 0.0;
   (void)context;

   for (i = 0; i < CMFEM_Vector_Size(x); i++)
   {
      const double value = CMFEM_Vector_Get(x, i);
      norm2 += value * value;
   }
   return exp(-30.0 * norm2);
}

static double initial_rate(const CMFEM_Vector *x, void *context)
{
   (void)x;
   (void)context;
   return 0.0;
}

static void wave_operator_reset_implicit(WaveOperatorContext *oper)
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

static void wave_operator_init(WaveOperatorContext *oper,
                               CMFEM_FiniteElementSpace *fespace,
                               const CMFEM_ArrayInt *ess_bdr,
                               double speed)
{
   CMFEM_ConstantCoefficient *one = NULL;

   oper->ess_tdof_list = CMFEM_ArrayInt_Construct();
   oper->mmat = CMFEM_SparseMatrix_Construct();
   oper->kmat = CMFEM_SparseMatrix_Construct();
   oper->z = CMFEM_Vector_ConstructSize(
                CMFEM_FiniteElementSpace_GetTrueVSize(fespace));
   oper->m_form = CMFEM_BilinearForm_New(fespace);
   oper->k_form = CMFEM_BilinearForm_New(fespace);
   oper->tmat = NULL;
   oper->m_prec = NULL;
   oper->t_prec = NULL;
   oper->c2 = CMFEM_ConstantCoefficient_New(speed * speed);

   CMFEM_BilinearForm_AddDomainIntegratorDiCc(oper->k_form, oper->c2);
   CMFEM_BilinearForm_Assemble(oper->k_form);

   one = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_BilinearForm_AddDomainIntegratorMiCc(oper->m_form, one);
   CMFEM_BilinearForm_Assemble(oper->m_form);
   CMFEM_ConstantCoefficient_Delete(one);

   CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace,
                                                 ess_bdr,
                                                 &oper->ess_tdof_list);
   CMFEM_BilinearForm_FormSystemMatrixSm(oper->k_form,
                                         &oper->ess_tdof_list,
                                         &oper->kmat);
   CMFEM_BilinearForm_FormSystemMatrixSm(oper->m_form,
                                         &oper->ess_tdof_list,
                                         &oper->mmat);
   oper->m_prec = CMFEM_DSmoother_NewSm(&oper->mmat);
}

static void wave_operator_destroy(WaveOperatorContext *oper)
{
   wave_operator_reset_implicit(oper);
   if (oper->m_prec != NULL)
   {
      CMFEM_DSmoother_Delete(oper->m_prec);
   }
   CMFEM_SparseMatrix_Destroy(&oper->kmat);
   CMFEM_SparseMatrix_Destroy(&oper->mmat);
   CMFEM_Vector_Destroy(&oper->z);
   CMFEM_ArrayInt_Destroy(&oper->ess_tdof_list);
   CMFEM_BilinearForm_Delete(oper->m_form);
   CMFEM_BilinearForm_Delete(oper->k_form);
   CMFEM_ConstantCoefficient_Delete(oper->c2);
}

static void wave_operator_mult(const CMFEM_Vector *u,
                               const CMFEM_Vector *du_dt,
                               CMFEM_Vector *d2udt2,
                               double time,
                               void *context)
{
   WaveOperatorContext *oper = (WaveOperatorContext *)context;
   (void)du_dt;
   (void)time;

   // Compute d2u_dt2 = M^{-1} (-K u).
   CMFEM_BilinearForm_FullMult(oper->k_form, u, &oper->z);
   CMFEM_Vector_Neg(&oper->z);
   CMFEM_Vector_SetSubVectorAi(&oper->z, &oper->ess_tdof_list, 0.0);
   CMFEM_PCGSmDs(&oper->mmat,
                 oper->m_prec,
                 &oper->z,
                 d2udt2,
                 0,
                 30,
                 1.0e-8,
                 0.0);
   CMFEM_Vector_SetSubVectorAi(d2udt2, &oper->ess_tdof_list, 0.0);
}

static void wave_operator_implicit_solve(double fac0,
                                         double fac1,
                                         const CMFEM_Vector *u,
                                         const CMFEM_Vector *du_dt,
                                         CMFEM_Vector *d2udt2,
                                         double time,
                                         void *context)
{
   WaveOperatorContext *oper = (WaveOperatorContext *)context;
   (void)fac1;
   (void)du_dt;
   (void)time;

   // Solve d2u_dt2 = M^{-1}[-K(u + fac0 * d2u_dt2)].
   if (oper->tmat == NULL)
   {
      oper->tmat = CMFEM_AddSmSm(1.0, &oper->mmat, fac0, &oper->kmat);
      oper->t_prec = CMFEM_DSmoother_NewSm(oper->tmat);
   }

   CMFEM_BilinearForm_FullMult(oper->k_form, u, &oper->z);
   CMFEM_Vector_Neg(&oper->z);
   CMFEM_Vector_SetSubVectorAi(&oper->z, &oper->ess_tdof_list, 0.0);
   CMFEM_PCGSmDs(oper->tmat,
                 oper->t_prec,
                 &oper->z,
                 d2udt2,
                 0,
                 100,
                 1.0e-8,
                 0.0);
   CMFEM_Vector_SetSubVectorAi(d2udt2, &oper->ess_tdof_list, 0.0);
}

static void wave_operator_set_parameters(WaveOperatorContext *oper,
                                         const CMFEM_Vector *u)
{
   (void)u;
   wave_operator_reset_implicit(oper);
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int ref_levels = 2;
   int order = 2;
   int ode_solver_type = 10;
   double t_final = 0.5;
   double dt = 1.0e-2;
   double speed = 1.0;
   int visualization = 1;
   int visit = 1;
   int dirichlet = 1;
   int vis_steps = 5;
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

      parsed = cmfem_parse_double_option(argc, argv, &i, "-c", "--speed",
                                         &speed);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-dir", "--dirichlet",
                                       "-neu", "--neumann",
                                       &dirichlet);
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
              "VisIt output is currently disabled in examples/c/ex23.c to avoid instability in the C port.\n");
      visit = 0;
   }

   // 2. Read the mesh and uniformly refine it.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   const int dim = CMFEM_Mesh_Dimension(mesh);

   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }

   // 3. Define the finite element space and the initial conditions.
   CMFEM_H1FeCollection *fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh,
                                                                          fec);
   CMFEM_GridFunction *u_gf = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction *dudt_gf = CMFEM_GridFunction_New(fespace);
   CMFEM_FunctionCoefficient *u0 =
      CMFEM_FunctionCoefficient_New(initial_solution, NULL);
   CMFEM_FunctionCoefficient *dudt0 =
      CMFEM_FunctionCoefficient_New(initial_rate, NULL);
   _Alignas(max_align_t) CMFEM_Vector u = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector dudt = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr = CMFEM_ArrayInt_Construct();
   WaveOperatorContext oper;
   CMFEM_SecondOrderTimeDependentOperator *oper_adapter = NULL;
   CMFEM_SecondOrderODESolver *ode_solver = NULL;
   CMFEM_VisItDataCollection *visit_dc = NULL;

   printf("Number of temperature unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   CMFEM_GridFunction_ProjectCoefficientFc(u_gf, u0);
   CMFEM_GridFunction_GetTrueDofs(u_gf, &u);

   CMFEM_GridFunction_ProjectCoefficientFc(dudt_gf, dudt0);
   CMFEM_GridFunction_GetTrueDofs(dudt_gf, &dudt);

   // 4. Mark the essential boundary attributes and initialize the wave
   //    operator context.
   if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
   {
      ess_bdr = CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(
                                                mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr, dirichlet ? 1 : 0);
   }

   wave_operator_init(&oper, fespace, &ess_bdr, speed);
   oper_adapter = CMFEM_SecondOrderTimeDependentOperator_New(
                     CMFEM_FiniteElementSpace_GetTrueVSize(fespace),
                     wave_operator_mult,
                     &oper,
                     wave_operator_implicit_solve,
                     &oper);
   ode_solver = CMFEM_SecondOrderODESolver_NewType(ode_solver_type);
   if (ode_solver == NULL)
   {
      fprintf(stderr, "Unsupported second-order ODE solver type: %d\n",
              ode_solver_type);
      CMFEM_SecondOrderTimeDependentOperator_Delete(oper_adapter);
      wave_operator_destroy(&oper);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_Vector_Destroy(&dudt);
      CMFEM_Vector_Destroy(&u);
      CMFEM_FunctionCoefficient_Delete(dudt0);
      CMFEM_FunctionCoefficient_Delete(u0);
      CMFEM_GridFunction_Delete(dudt_gf);
      CMFEM_GridFunction_Delete(u_gf);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_H1FeCollection_Delete(fec);
      CMFEM_Mesh_Delete(mesh);
      return 1;
   }

   // 5. Save the initial state and initialize output.
   CMFEM_GridFunction_SetFromTrueDofs(u_gf, &u);
   CMFEM_GridFunction_SetFromTrueDofs(dudt_gf, &dudt);
   CMFEM_Mesh_Print(mesh, "ex23.mesh", precision);
   CMFEM_GridFunction_Save(u_gf, "ex23-init.gf", precision);
   CMFEM_GridFunction_Save(dudt_gf, "ex23-init-rate.gf", precision);

   if (visit)
   {
      visit_dc = CMFEM_VisItDataCollection_New("Example23", mesh);
      CMFEM_VisItDataCollection_RegisterFieldGf(visit_dc, "solution", u_gf);
      CMFEM_VisItDataCollection_RegisterFieldGf(visit_dc, "rate", dudt_gf);
      CMFEM_VisItDataCollection_SetCycle(visit_dc, 0);
      CMFEM_VisItDataCollection_SetTime(visit_dc, 0.0);
      CMFEM_VisItDataCollection_Save(visit_dc);
   }

   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, u_gf, "localhost", 19916);
   }

   // 6. Perform time integration with the selected second-order ODE solver.
   {
      double t = 0.0;
      int ti = 0;
      int last_step = 0;

      CMFEM_SecondOrderODESolver_Init(ode_solver, oper_adapter);

      while (!last_step)
      {
         ti++;
         if (t + dt >= t_final - dt / 2.0)
         {
            last_step = 1;
         }

         CMFEM_SecondOrderODESolver_Step(ode_solver, &u, &dudt, &t, &dt);

         if (last_step || (ti % vis_steps) == 0)
         {
            printf("step %d, t = %.17g\n", ti, t);

            CMFEM_GridFunction_SetFromTrueDofs(u_gf, &u);
            CMFEM_GridFunction_SetFromTrueDofs(dudt_gf, &dudt);

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

         wave_operator_set_parameters(&oper, &u);
      }
   }

   // 7. Save the final solution for later visualization.
   CMFEM_GridFunction_SetFromTrueDofs(u_gf, &u);
   CMFEM_GridFunction_SetFromTrueDofs(dudt_gf, &dudt);
   CMFEM_GridFunction_Save(u_gf, "ex23-final.gf", precision);
   CMFEM_GridFunction_Save(dudt_gf, "ex23-final-rate.gf", precision);

   // 8. Free the used memory.
   if (visit_dc != NULL)
   {
      CMFEM_VisItDataCollection_Delete(visit_dc);
   }
   CMFEM_SecondOrderODESolver_Delete(ode_solver);
   CMFEM_SecondOrderTimeDependentOperator_Delete(oper_adapter);
   wave_operator_destroy(&oper);
   CMFEM_ArrayInt_Destroy(&ess_bdr);
   CMFEM_Vector_Destroy(&dudt);
   CMFEM_Vector_Destroy(&u);
   CMFEM_FunctionCoefficient_Delete(dudt0);
   CMFEM_FunctionCoefficient_Delete(u0);
   CMFEM_GridFunction_Delete(dudt_gf);
   CMFEM_GridFunction_Delete(u_gf);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_H1FeCollection_Delete(fec);
   CMFEM_Mesh_Delete(mesh);

   return 0;
}
