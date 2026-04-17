//                              MFEM C Example 3
//
// Description: This example code solves a simple electromagnetic diffusion
//              problem corresponding to the second order definite Maxwell
//              equation curl curl E + E = f with boundary condition E x n
//              equal to a given tangential field.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

static double cmfem_ex3_freq = 1.0;
static double cmfem_ex3_kappa = 0.0;
static int cmfem_ex3_dim = 0;
static const double cmfem_ex3_pi = 3.14159265358979323846;

// Exact solution, E, and right-hand side, f.
static void cmfem_ex3_E_exact(const CMFEM_Vector *x, CMFEM_Vector *E, void *context)
{
   (void)context;
   if (cmfem_ex3_dim == 3)
   {
      CMFEM_Vector_Set(E, 0, sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 1)));
      CMFEM_Vector_Set(E, 1, sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 2)));
      CMFEM_Vector_Set(E, 2, sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 0)));
   }
   else
   {
      CMFEM_Vector_Set(E, 0, sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 1)));
      CMFEM_Vector_Set(E, 1, sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 0)));
      if (CMFEM_Vector_Size(x) == 3)
      {
         CMFEM_Vector_Set(E, 2, 0.0);
      }
   }
}

static void cmfem_ex3_f_exact(const CMFEM_Vector *x, CMFEM_Vector *f, void *context)
{
   double factor = 1.0 + cmfem_ex3_kappa * cmfem_ex3_kappa;
   (void)context;
   if (cmfem_ex3_dim == 3)
   {
      CMFEM_Vector_Set(f, 0, factor * sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 1)));
      CMFEM_Vector_Set(f, 1, factor * sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 2)));
      CMFEM_Vector_Set(f, 2, factor * sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 0)));
   }
   else
   {
      CMFEM_Vector_Set(f, 0, factor * sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 1)));
      CMFEM_Vector_Set(f, 1, factor * sin(cmfem_ex3_kappa * CMFEM_Vector_Get(x, 0)));
      if (CMFEM_Vector_Size(x) == 3)
      {
         CMFEM_Vector_Set(f, 2, 0.0);
      }
   }
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("beam-tet.mesh");
   int order = 1;
   int static_cond = 0;
   int pa = 0;
   int nc = 0;
   const char *device_config = "cpu";
   int visualization = 1;
   int sdim;
   int ref_levels;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh", &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-f", "--frequency", &cmfem_ex3_freq);
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

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-nc", "--non-conforming",
                                       "-c", "--conforming",
                                       &nc);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device", &device_config);
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

   cmfem_ex3_kappa = cmfem_ex3_freq * cmfem_ex3_pi;

   // 2. Enable hardware devices and programming models based on command line
   //    options.
   CMFEM_Device *device = CMFEM_Device_New(device_config);
   CMFEM_Device_Print(device);

   // 3. Read the mesh from the given mesh file.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   cmfem_ex3_dim = CMFEM_Mesh_Dimension(mesh);
   sdim = CMFEM_Mesh_SpaceDimension(mesh);
   if (nc)
   {
      CMFEM_Mesh_EnsureNCMesh(mesh, 1);
   }

   // 4. Refine the mesh to increase the resolution.
   ref_levels = cmfem_uniform_refinement_levels(50000.0, CMFEM_Mesh_GetNE(mesh), cmfem_ex3_dim);
   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }

   // 5. Define a Nedelec finite element space on the mesh.
   CMFEM_ND_FECollection *fec = CMFEM_ND_FECollection_NewOrderDim(order, cmfem_ex3_dim);
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshND(mesh, fec);
   printf("Number of finite element unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   // 6. Determine the list of true essential boundary dofs.
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
   if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
   {
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
         CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr, 1);
      CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr, &ess_tdof_list);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side
   //    of the discrete system.
   CMFEM_VectorFunctionCoefficient *f =
      CMFEM_VectorFunctionCoefficient_New(sdim, cmfem_ex3_f_exact, NULL);
   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_LinearForm_AddDomainIntegrator_VectorFEDomainLFIntegrator(b, f);
   CMFEM_LinearForm_Assemble(b);

   // 8. Define the solution vector x as a finite element grid function and
   //    initialize it by projecting the exact solution.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_VectorFunctionCoefficient *E =
      CMFEM_VectorFunctionCoefficient_New(sdim, cmfem_ex3_E_exact, NULL);
   CMFEM_GridFunction_ProjectVectorFunctionCoefficient(x, E);

   // 9. Set up the bilinear form corresponding to the EM diffusion operator.
   CMFEM_ConstantCoefficient *muinv = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_ConstantCoefficient *sigma = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   if (pa)
   {
      CMFEM_BilinearForm_SetAssemblyLevelPartial(a);
   }
   CMFEM_BilinearForm_AddDomainIntegrator_CurlCurl(a, muinv);
   CMFEM_BilinearForm_AddDomainIntegrator_VectorFEMass(a, sigma);
   if (static_cond)
   {
      CMFEM_BilinearForm_EnableStaticCondensation(a);
   }
   CMFEM_BilinearForm_Assemble(a);

   // 10. Assemble the linear system A X = B.
   _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
   CMFEM_BilinearForm_FormLinearSystemOperator(a, &ess_tdof_list, x, b, &A, &X, &B);

   printf("Size of linear system: %d\n", CMFEM_OperatorPtr_Height(&A));

   // 11. Solve the linear system A X = B.
   if (pa)
   {
      CMFEM_OperatorJacobiSmoother *M =
         CMFEM_OperatorJacobiSmoother_NewBilinearForm(a, &ess_tdof_list);
      CMFEM_PCG_OperatorJacobiSmoother(&A, M, &B, &X, 1, 1000, 1e-12, 0.0);
      CMFEM_OperatorJacobiSmoother_Delete(M);
   }
   else
   {
      CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewOperator(&A);
      CMFEM_PCG_OperatorGSSmoother(&A, M, &B, &X, 1, 500, 1e-12, 0.0);
      CMFEM_GSSmoother_Delete(M);
   }

   // 12. Recover the solution as a finite element grid function.
   CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
   // 13. Compute and print the L^2 norm of the error.
   printf("\n|| E_h - E ||_{L^2} = %.16g\n\n",
          CMFEM_GridFunction_ComputeL2ErrorVectorFunctionCoefficient(x, E));

   // 14. Save the refined mesh and the solution.
   CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
   CMFEM_GridFunction_Save(x, "sol.gf", 8);

   // 15. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
   }

   // 16. Free the used memory.
   CMFEM_OperatorPtr_Destroy(&A);
   CMFEM_Vector_Destroy(&B);
   CMFEM_Vector_Destroy(&X);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_ConstantCoefficient_Delete(sigma);
   CMFEM_ConstantCoefficient_Delete(muinv);
   CMFEM_GridFunction_Delete(x);
   CMFEM_VectorFunctionCoefficient_Delete(E);
   CMFEM_LinearForm_Delete(b);
   CMFEM_VectorFunctionCoefficient_Delete(f);
   CMFEM_ArrayInt_Destroy(&ess_tdof_list);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_ND_FECollection_Delete(fec);
   CMFEM_Mesh_Delete(mesh);
   CMFEM_Device_Delete(device);
   return 0;
}
