//                              MFEM C Example 18
//
// Description: This example solves the compressible Euler equations with a DG
//              spatial discretization and explicit time integration.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>
#include <time.h>

struct ex18_context
{
   void *dghcl;
};

static void ex18_mult_t(const CMFEM_Vector *x,
                        CMFEM_Vector *y,
                        double time,
                        void *context)
{
   struct ex18_context *ctx = (struct ex18_context *)context;
   (void)time;
   CMFEM_Ex18_MultDghcl(ctx->dghcl, x, y);
}

static double elapsed_seconds(clock_t start, clock_t end)
{
   return (double)(end - start) / (double)CLOCKS_PER_SEC;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int problem = 1;
   const double specific_heat_ratio = 1.4;
   const double gas_constant = 1.0;
   const char *mesh_file = NULL;
   int int_order_offset = 1;
   int ref_levels = 1;
   int order = 3;
   int ode_solver_type = 4;
   double t_final = 2.0;
   double dt = -0.01;
   double cfl = 0.3;
   int visualization = 1;
   int preassemble_weak_div = 1;
   int vis_steps = 50;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-p", "--problem",
                                      &problem);
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

      parsed = cmfem_parse_double_option(argc, argv, &i, "-c", "--cfl-number",
                                         &cfl);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-ea", "--element-assembly-divergence",
                                       "-mf", "--matrix-free-divergence",
                                       &preassemble_weak_div);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-vs",
                                      "--visualization-steps", &vis_steps);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   if (problem < 1 || problem > 4)
   {
      fprintf(stderr, "Unrecognized problem type: %d\n", problem);
      return 1;
   }

   // 2. Read and refine the mesh.
   {
      CMFEM_Mesh *mesh = mesh_file == NULL ?
                         CMFEM_Ex18_NewDefaultMesh(problem) :
                         CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      const int dim = CMFEM_Mesh_Dimension(mesh);
      const int num_equations = dim + 2;
      CMFEM_ODESolver *ode_solver = NULL;
      CMFEM_DgFeCollection *fec = NULL;
      CMFEM_FiniteElementSpace *fes = NULL;
      CMFEM_FiniteElementSpace *dfes = NULL;
      CMFEM_FiniteElementSpace *vfes = NULL;
      CMFEM_GridFunction *sol = NULL;
      CMFEM_TimeDependentOperator *oper = NULL;
      _Alignas(max_align_t) CMFEM_Vector state = CMFEM_Vector_Construct();
      struct ex18_context ctx;
      clock_t start_time;
      clock_t end_time;
      double t = 0.0;
      int ti = 0;
      int done = 0;

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      // 3. Define the ODE solver and DG finite element spaces.
      ode_solver = CMFEM_ODESolver_NewType(ode_solver_type);
      fec = CMFEM_DgFeCollection_NewOrderDimBasis(order,
                                                  dim,
                                                  CMFEM_BASIS_GAUSS_LEGENDRE);
      fes = CMFEM_FiniteElementSpace_NewMeshDg(mesh, fec);
      dfes = CMFEM_FiniteElementSpace_NewMeshDgVDim(mesh, fec, dim);
      vfes = CMFEM_FiniteElementSpace_NewMeshDgVDim(mesh, fec, num_equations);

      printf("Number of unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetVSize(vfes));

      // 4. Project the initial condition and save the initial state.
      sol = CMFEM_GridFunction_New(vfes);
      CMFEM_Ex18_ProjectInitialConditionGf(problem,
                                           specific_heat_ratio,
                                           gas_constant,
                                           sol);
      CMFEM_Ex18_SaveStateMeshSfesGf(mesh,
                                     fes,
                                     num_equations,
                                     sol,
                                     "euler-mesh.mesh",
                                     "init");

      // 5. Set up the DG Euler operator and visualization.
      ctx.dghcl = CMFEM_Ex18_NewDghcl(vfes,
                                      int_order_offset,
                                      specific_heat_ratio,
                                      preassemble_weak_div);
      if (visualization)
      {
         if (!CMFEM_Ex18_SendMomentumToGLVisMeshSfesDfesGf(mesh,
                                                           fes,
                                                           dfes,
                                                           sol,
                                                           0.0,
                                                           1))
         {
            visualization = 0;
            printf("Unable to connect to GLVis server at localhost:19916\n");
            printf("GLVis visualization disabled.\n");
         }
         else
         {
            printf("GLVis visualization paused. Press space (in the GLVis window) to resume it.\n");
         }
      }

      if (cfl > 0.0)
      {
         dt = CMFEM_Ex18_ComputeStableDtMeshDghclGf(mesh,
                                                    order,
                                                    cfl,
                                                    ctx.dghcl,
                                                    sol);
      }

      CMFEM_GridFunction_CopyToVec(sol, &state);
      oper = CMFEM_TimeDependentOperator_New(CMFEM_Vector_Size(&state),
                                             CMFEM_TIME_DEPENDENT_OPERATOR_EXPLICIT,
                                             ex18_mult_t,
                                             &ctx,
                                             NULL,
                                             NULL);
      CMFEM_Ex18_SetTimeDghcl(ctx.dghcl, t);
      CMFEM_ODESolver_Init(ode_solver, oper);

      start_time = clock();

      // 6. Integrate in time.
      while (!done)
      {
         double dt_real = dt;
         if (dt_real > t_final - t)
         {
            dt_real = t_final - t;
         }

         CMFEM_ODESolver_Step(ode_solver, &state, &t, &dt_real);
         CMFEM_GridFunction_SetFromVec(sol, &state);

         if (cfl > 0.0)
         {
            dt = CMFEM_Ex18_ComputeStableDtMeshDghclGf(mesh,
                                                       order,
                                                       cfl,
                                                       ctx.dghcl,
                                                       sol);
         }

         ti++;
         done = (t >= t_final - 1.0e-8 * dt);
         if (done || ti % vis_steps == 0)
         {
            printf("time step: %d, time: %.8g\n", ti, t);
            if (visualization)
            {
               if (!CMFEM_Ex18_SendMomentumToGLVisMeshSfesDfesGf(mesh,
                                                                 fes,
                                                                 dfes,
                                                                 sol,
                                                                 t,
                                                                 0))
               {
                  visualization = 0;
                  printf("GLVis visualization disabled.\n");
               }
            }
         }
      }

      end_time = clock();
      printf(" done, %gs.\n", elapsed_seconds(start_time, end_time));

      // 7. Save the final state and report the final error.
      CMFEM_Ex18_SaveStateMeshSfesGf(mesh,
                                     fes,
                                     num_equations,
                                     sol,
                                     "euler-mesh-final.mesh",
                                     "final");
      printf("Solution error: %.8g\n",
             CMFEM_Ex18_ComputeL2ErrorGf(problem,
                                         specific_heat_ratio,
                                         gas_constant,
                                         sol));

      CMFEM_TimeDependentOperator_Delete(oper);
      CMFEM_Ex18_DeleteDghcl(ctx.dghcl);
      CMFEM_Vector_Destroy(&state);
      CMFEM_GridFunction_Delete(sol);
      CMFEM_FiniteElementSpace_Delete(vfes);
      CMFEM_FiniteElementSpace_Delete(dfes);
      CMFEM_FiniteElementSpace_Delete(fes);
      CMFEM_DgFeCollection_Delete(fec);
      CMFEM_ODESolver_Delete(ode_solver);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
