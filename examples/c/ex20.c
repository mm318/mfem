//                              MFEM C Example 20
//
// Description: This example demonstrates MFEM's symplectic time integrator on
//              a small family of one-dimensional Hamiltonian systems.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int prob_ = 0;
static double m_ = 1.0;
static double k_ = 1.0;

static double hamiltonian(double q, double p, double t)
{
   double h = 1.0 - 0.5 / m_ + 0.5 * p * p / m_;
   (void)t;

   switch (prob_)
   {
      case 1:
         h += k_ * (1.0 - cos(q));
         break;
      case 2:
         h += k_ * (1.0 - exp(-0.5 * q * q));
         break;
      case 3:
         h += 0.5 * k_ * (1.0 + q * q) * q * q;
         break;
      case 4:
         h += 0.5 * k_ * (1.0 - 0.125 * q * q) * q * q;
         break;
      default:
         h += 0.5 * k_ * q * q;
         break;
   }
   return h;
}

static void grad_t(const CMFEM_Vector *x, CMFEM_Vector *y, void *context)
{
   (void)context;
   CMFEM_Vector_Set(y, 0, CMFEM_Vector_Get(x, 0) / m_);
}

static void neg_grad_v(const CMFEM_Vector *x,
                       CMFEM_Vector *y,
                       double time,
                       void *context)
{
   const double q = CMFEM_Vector_Get(x, 0);
   (void)time;
   (void)context;

   switch (prob_)
   {
      case 1:
         CMFEM_Vector_Set(y, 0, -k_ * sin(q));
         break;
      case 2:
         CMFEM_Vector_Set(y, 0, -k_ * q * exp(-0.5 * q * q));
         break;
      case 3:
         CMFEM_Vector_Set(y, 0, -k_ * (1.0 + 2.0 * q * q) * q);
         break;
      case 4:
         CMFEM_Vector_Set(y, 0, -k_ * (1.0 - 0.25 * q * q) * q);
         break;
      default:
         CMFEM_Vector_Set(y, 0, -k_ * q);
         break;
   }
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int order = 1;
   int nsteps = 100;
   double dt = 0.1;
   int visualization = 1;
   int gnuplot = 0;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order",
                                          &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-p", "--problem-type",
                                      &prob_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-n", "--number-of-steps",
                                      &nsteps);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-dt", "--time-step",
                                         &dt);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-m", "--mass", &m_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-k", "--spring-const",
                                         &k_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-gp", "--gnuplot",
                                       "-no-gp", "--no-gnuplot",
                                       &gnuplot);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   if (nsteps < 1)
   {
      fprintf(stderr, "Number of time steps must be positive.\n");
      return 1;
   }

   // 2. Create and initialize the symplectic solver.
   CMFEM_SIAVSolver *solver = CMFEM_SIAVSolver_NewOrder(order);
   CMFEM_SIAVSolver_Init(solver, 1, grad_t, NULL, neg_grad_v, NULL);

   // 3. Set the initial conditions for the canonical state.
   double t = 0.0;
   CMFEM_Vector *q = CMFEM_Vector_NewSize(1);
   CMFEM_Vector *p = CMFEM_Vector_NewSize(1);
   double *e = (double *)malloc((size_t)(nsteps + 1) * sizeof(double));
   FILE *ofs = NULL;

   if (e == NULL)
   {
      fprintf(stderr, "Unable to allocate the energy history.\n");
      CMFEM_Vector_Delete(p);
      CMFEM_Vector_Delete(q);
      CMFEM_SIAVSolver_Delete(solver);
      return 1;
   }

   CMFEM_Vector_Set(q, 0, 0.0);
   CMFEM_Vector_Set(p, 0, 1.0);

   // 4. Prepare GnuPlot output if requested.
   if (gnuplot)
   {
      ofs = fopen("ex20.dat", "w");
      if (ofs == NULL)
      {
         fprintf(stderr, "Unable to open ex20.dat for writing.\n");
         free(e);
         CMFEM_Vector_Delete(p);
         CMFEM_Vector_Delete(q);
         CMFEM_SIAVSolver_Delete(solver);
         return 1;
      }
      fprintf(ofs, "%.17g\t%.17g\t%.17g\n", t,
              CMFEM_Vector_Get(q, 0), CMFEM_Vector_Get(p, 0));
   }

   // 5. Create the ribbon mesh used to visualize the trajectory in phase
   //    space.
   CMFEM_Mesh *mesh = NULL;
   if (visualization)
   {
      mesh = CMFEM_Mesh_NewDimensionVerticesElementsBoundarySpace(
                2, 2 * (nsteps + 1), nsteps, 0, 3);
   }

   // 6. Advance the Hamiltonian system and record its state.
   {
      double e_mean = 0.0;

      for (i = 0; i < nsteps; i++)
      {
         if (i == 0)
         {
            double x0[3] = { 0.0, 0.0, 0.0 };
            double x1[3] = { CMFEM_Vector_Get(q, 0), CMFEM_Vector_Get(p, 0), 0.0 };
            e[0] = hamiltonian(CMFEM_Vector_Get(q, 0), CMFEM_Vector_Get(p, 0), t);
            e_mean += e[0];

            if (visualization)
            {
               CMFEM_Mesh_AddVertex(mesh, x0);
               CMFEM_Mesh_AddVertex(mesh, x1);
            }
         }

         CMFEM_SIAVSolver_Step(solver, q, p, &t, &dt);
         e[i + 1] = hamiltonian(CMFEM_Vector_Get(q, 0), CMFEM_Vector_Get(p, 0), t);
         e_mean += e[i + 1];

         if (gnuplot)
         {
            fprintf(ofs, "%.17g\t%.17g\t%.17g\t%.17g\n",
                    t,
                    CMFEM_Vector_Get(q, 0),
                    CMFEM_Vector_Get(p, 0),
                    e[i + 1]);
         }

         if (visualization)
         {
            double x0[3] = { 0.0, 0.0, t };
            double x1[3] = { CMFEM_Vector_Get(q, 0), CMFEM_Vector_Get(p, 0), t };
            int v[4] = { 2 * i, 2 * (i + 1), 2 * (i + 1) + 1, 2 * i + 1 };
            CMFEM_Mesh_AddVertex(mesh, x0);
            CMFEM_Mesh_AddVertex(mesh, x1);
            CMFEM_Mesh_AddQuad(mesh, v, 1);
         }
      }

      // 7. Compute and display the mean and standard deviation of the energy.
      {
         double e_var = 0.0;
         double e_sd;
         e_mean /= (double)(nsteps + 1);
         for (i = 0; i <= nsteps; i++)
         {
            e_var += (e[i] - e_mean) * (e[i] - e_mean);
         }
         e_var /= (double)(nsteps + 1);
         e_sd = sqrt(e_var);
         printf("\nMean and standard deviation of the energy\n");
         printf("%.17g\t%.17g\n", e_mean, e_sd);
      }
   }

   // 8. Finalize the GnuPlot output files.
   if (gnuplot)
   {
      fclose(ofs);
      ofs = fopen("gnuplot_ex20.inp", "w");
      if (ofs != NULL)
      {
         fprintf(ofs,
                 "plot 'ex20.dat' using 1:2 w l t 'q', 'ex20.dat' using 1:3 w l t 'p', 'ex20.dat' using 1:4 w l t 'H'\n");
         fclose(ofs);
      }
   }

   // 9. Finalize the phase-space visualization.
   if (visualization)
   {
      CMFEM_H1FeCollection *fec;
      CMFEM_FiniteElementSpace *fespace;
      CMFEM_GridFunction *energy;

      CMFEM_Mesh_FinalizeQuadMesh(mesh, 1, 0, 0);
      fec = CMFEM_H1FeCollection_NewOrderDim(1, 2);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, fec);
      energy = CMFEM_GridFunction_New(fespace);
      CMFEM_GridFunction_Assign(energy, 0.0);

      for (i = 0; i <= nsteps; i++)
      {
         CMFEM_GridFunction_Set(energy, 2 * i + 0, e[i]);
         CMFEM_GridFunction_Set(energy, 2 * i + 1, e[i]);
      }

      CMFEM_SendSolutionToGLVis(mesh, energy, "localhost", 19916);
      CMFEM_GridFunction_Delete(energy);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_H1FeCollection_Delete(fec);
      CMFEM_Mesh_Delete(mesh);
   }

   // 10. Free the remaining resources.
   free(e);
   CMFEM_Vector_Delete(p);
   CMFEM_Vector_Delete(q);
   CMFEM_SIAVSolver_Delete(solver);
   return 0;
}
