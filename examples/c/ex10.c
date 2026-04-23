//                              MFEM C Example 10
//
// Description: This example solves a time-dependent nonlinear elasticity
//              problem using a focused CMFEM bridge for the hyperelastic
//              operator and MFEM's ODE solvers through the C wrapper layer.

#include "common.h"
#include "cmfem.h"
#include "adapters/ex10/hyperelastic_dynamics.h"

#include <math.h>
#include <stdio.h>

typedef struct Ex10Context
{
   void *oper;
} Ex10Context;

static void initial_deformation(const CMFEM_Vector *x,
                                CMFEM_Vector *y,
                                void *context)
{
   int i;
   (void)context;
   for (i = 0; i < CMFEM_Vector_Size(x); i++)
   {
      CMFEM_Vector_Set(y, i, CMFEM_Vector_Get(x, i));
   }
}

static void initial_velocity(const CMFEM_Vector *x,
                             CMFEM_Vector *v,
                             void *context)
{
   const int dim = CMFEM_Vector_Size(x);
   const double s = 0.1 / 64.0;
   int i;
   (void)context;

   for (i = 0; i < dim; i++)
   {
      CMFEM_Vector_Set(v, i, 0.0);
   }
   CMFEM_Vector_Set(v,
                    dim - 1,
                    s * CMFEM_Vector_Get(x, 0) * CMFEM_Vector_Get(x, 0) *
                    (8.0 - CMFEM_Vector_Get(x, 0)));
   CMFEM_Vector_Set(v,
                    0,
                    -s * CMFEM_Vector_Get(x, 0) * CMFEM_Vector_Get(x, 0));
}

static void ex10_mult(const CMFEM_Vector *state,
                      CMFEM_Vector *dstate_dt,
                      double time,
                      void *context)
{
   Ex10Context *ctx = (Ex10Context *)context;
   (void)time;
   CMFEM_Ex10HyperelasticOperator_Mult(ctx->oper, state, dstate_dt);
}

static void ex10_implicit_solve(double gamma,
                                const CMFEM_Vector *state,
                                CMFEM_Vector *k,
                                double time,
                                void *context)
{
   Ex10Context *ctx = (Ex10Context *)context;
   (void)time;
   CMFEM_Ex10HyperelasticOperator_ImplicitSolve(ctx->oper, gamma, state, k);
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("beam-quad.mesh");
   int ref_levels = 2;
   int order = 2;
   int ode_solver_type = 23;
   double t_final = 300.0;
   double dt = 3.0;
   double viscosity = 1.0e-2;
   double mu = 0.25;
   double bulk_modulus = 5.0;
   int visualization = 1;
   int vis_steps = 1;
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

      parsed = cmfem_parse_double_option(argc, argv, &i, "-v", "--viscosity",
                                         &viscosity);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-mu",
                                         "--shear-modulus", &mu);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-K",
                                         "--bulk-modulus", &bulk_modulus);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-vs",
                                      "--visualization-steps", &vis_steps);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Read the mesh and define the requested ODE solver.
   {
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      const int dim = CMFEM_Mesh_Dimension(mesh);
      CMFEM_ODESolver *ode_solver = CMFEM_ODESolver_NewType(ode_solver_type);
      CMFEM_H1FeCollection *fe_coll = NULL;
      CMFEM_FiniteElementSpace *fespace = NULL;
      CMFEM_GridFunction *v = NULL;
      CMFEM_GridFunction *x = NULL;
      CMFEM_GridFunction *x_ref = NULL;
      CMFEM_GridFunction *x_def = NULL;
      CMFEM_L2FeCollection *w_fec = NULL;
      CMFEM_FiniteElementSpace *w_fespace = NULL;
      CMFEM_GridFunction *w = NULL;
      CMFEM_VectorFunctionCoefficient *velo = NULL;
      CMFEM_VectorFunctionCoefficient *deform = NULL;
      Ex10Context ctx;
      CMFEM_TimeDependentOperator *oper_adapter = NULL;

      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_Vector state = CMFEM_Vector_Construct();

      if (ode_solver == NULL)
      {
         fprintf(stderr, "Unsupported ODE solver type: %d\n", ode_solver_type);
         CMFEM_Vector_Destroy(&state);
         CMFEM_ArrayInt_Destroy(&ess_bdr);
         CMFEM_Mesh_Delete(mesh);
         return 1;
      }

      // 3. Refine the mesh to increase the resolution.
      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      // 4. Define the vector finite element space for the deformation and
      //    velocity, as well as the discontinuous space for the energy
      //    density.
      fe_coll = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1VDim(mesh, fe_coll, dim);
      v = CMFEM_GridFunction_New(fespace);
      x = CMFEM_GridFunction_New(fespace);
      x_ref = CMFEM_GridFunction_New(fespace);
      x_def = CMFEM_GridFunction_New(fespace);
      CMFEM_Mesh_GetNodesGf(mesh, x_ref);

      w_fec = CMFEM_L2FeCollection_NewOrderDim(order + 1, dim);
      w_fespace = CMFEM_FiniteElementSpace_NewMeshL2Vdim(mesh, w_fec, 1);
      w = CMFEM_GridFunction_New(w_fespace);

      printf("Number of velocity/deformation unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

      // 5. Set the initial conditions and the fixed essential boundary.
      velo = CMFEM_VectorFunctionCoefficient_New(dim, initial_velocity, NULL);
      deform = CMFEM_VectorFunctionCoefficient_New(dim, initial_deformation, NULL);
      CMFEM_GridFunction_ProjectCoefficientVfc(v, velo);
      CMFEM_GridFunction_SetTrueVector(v);
      CMFEM_GridFunction_ProjectCoefficientVfc(x, deform);
      CMFEM_GridFunction_SetTrueVector(x);

      ess_bdr = CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(
                                                mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr, 0);
      if (CMFEM_ArrayInt_Size(&ess_bdr) > 0)
      {
         CMFEM_ArrayInt_Set(&ess_bdr, 0, 1);
      }

      // 6. Initialize the hyperelastic operator, pack the initial state, and
      //    print the initial energies.
      ctx.oper = CMFEM_Ex10HyperelasticOperator_New(fespace,
                                                    &ess_bdr,
                                                    viscosity,
                                                    mu,
                                                    bulk_modulus);
      CMFEM_Ex10HyperelasticOperator_SetState(ctx.oper, v, x, &state);
      printf("initial elastic energy (EE) = %.17g\n",
             CMFEM_Ex10HyperelasticOperator_ElasticEnergy(ctx.oper, &state));
      printf("initial kinetic energy (KE) = %.17g\n",
             CMFEM_Ex10HyperelasticOperator_KineticEnergy(ctx.oper, &state));
      printf("initial   total energy (TE) = %.17g\n",
             CMFEM_Ex10HyperelasticOperator_ElasticEnergy(ctx.oper, &state) +
             CMFEM_Ex10HyperelasticOperator_KineticEnergy(ctx.oper, &state));

      oper_adapter = CMFEM_TimeDependentOperator_New(
                        CMFEM_Ex10HyperelasticOperator_StateSize(ctx.oper),
                        CMFEM_TIME_DEPENDENT_OPERATOR_IMPLICIT,
                        ex10_mult,
                        &ctx,
                        ex10_implicit_solve,
                        &ctx);
      CMFEM_TimeDependentOperator_SetImplicitVariableType(
         oper_adapter,
         CMFEM_IMPLICIT_VARIABLE_SLOPE);
      CMFEM_ODESolver_Init(ode_solver, oper_adapter);

      // 7. Perform time integration, reporting the total energy drift.
      {
         double t = 0.0;
         const double initial_total_energy =
            CMFEM_Ex10HyperelasticOperator_ElasticEnergy(ctx.oper, &state) +
            CMFEM_Ex10HyperelasticOperator_KineticEnergy(ctx.oper, &state);
         int ti;
         int last_step = 0;

         if (visualization)
         {
            printf("ex10.c uses a simplified undeformed-mesh GLVis path.\n");
            CMFEM_Ex10HyperelasticOperator_SetFields(ctx.oper, &state, v, x);
            CMFEM_SendSolutionToGLVis(mesh, v, "localhost", 19916);
         }

         for (ti = 1; !last_step; ti++)
         {
            double dt_real = fmin(dt, t_final - t);
            double ee;
            double ke;

            CMFEM_ODESolver_Step(ode_solver, &state, &t, &dt_real);
            last_step = (t >= t_final - 1.0e-8 * dt);

            if (last_step || (ti % vis_steps) == 0)
            {
               ee = CMFEM_Ex10HyperelasticOperator_ElasticEnergy(ctx.oper,
                                                                 &state);
               ke = CMFEM_Ex10HyperelasticOperator_KineticEnergy(ctx.oper,
                                                                 &state);

               printf("step %d, t = %.17g, EE = %.17g, KE = %.17g, dTE = %.17g\n",
                      ti,
                      t,
                      ee,
                      ke,
                      ee + ke - initial_total_energy);

               if (visualization)
               {
                  CMFEM_Ex10HyperelasticOperator_SetFields(ctx.oper, &state, v, x);
                  CMFEM_SendSolutionToGLVis(mesh, v, "localhost", 19916);
               }
            }
         }
      }

      // 8. Save the displaced mesh, the final velocity, and the elastic
      //    energy density.
      CMFEM_Ex10HyperelasticOperator_SetFields(ctx.oper, &state, v, x);
      CMFEM_Ex10HyperelasticOperator_GetElasticEnergyDensity(ctx.oper, &state, w);
      CMFEM_GridFunction_Save(v, "velocity.sol", 8);
      CMFEM_GridFunction_Save(w, "elastic_energy.sol", 8);

      CMFEM_GridFunction_CopyFromGf(x_def, x);
      CMFEM_GridFunction_AddScaledGf(x_def, -1.0, x_ref);
      CMFEM_Mesh_SetNodalFESpace(mesh, fespace);
      CMFEM_Mesh_AddDisplacementToNodes(mesh, x_def);
      CMFEM_Mesh_Print(mesh, "deformed.mesh", 8);

      // 9. Free the used memory.
      CMFEM_TimeDependentOperator_Delete(oper_adapter);
      CMFEM_Ex10HyperelasticOperator_Delete(ctx.oper);
      CMFEM_VectorFunctionCoefficient_Delete(deform);
      CMFEM_VectorFunctionCoefficient_Delete(velo);
      CMFEM_GridFunction_Delete(w);
      CMFEM_FiniteElementSpace_Delete(w_fespace);
      CMFEM_L2FeCollection_Delete(w_fec);
      CMFEM_GridFunction_Delete(x_def);
      CMFEM_GridFunction_Delete(x_ref);
      CMFEM_GridFunction_Delete(x);
      CMFEM_GridFunction_Delete(v);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_H1FeCollection_Delete(fe_coll);
      CMFEM_Vector_Destroy(&state);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_ODESolver_Delete(ode_solver);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
