//                              MFEM C Example 21
//
// Description: This example extends Example 2 with an adaptive refinement loop
//              for the linear elasticity cantilever beam problem.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

static CMFEM_Mesh *make_deformed_mesh_copy(CMFEM_Mesh *mesh,
                                           CMFEM_FiniteElementSpace *fespace,
                                           const CMFEM_GridFunction *x)
{
   CMFEM_Mesh *deformed;

   if (!CMFEM_Mesh_HasNodes(mesh))
   {
      CMFEM_Mesh_SetNodalFESpace(mesh, fespace);
   }

   deformed = CMFEM_Mesh_NewCopy(mesh);
   CMFEM_Mesh_AddDisplacementToNodes(deformed, x);
   return deformed;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("beam-tri.mesh");
   int order = 1;
   int static_cond = 0;
   int flux_averaging = 0;
   int visualization = 1;
   int dim;
   int tdim;
   int it;

   for (it = 1; it < argc; it++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &it, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &it, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&it, argv,
                                       "-sc", "--static-condensation",
                                       "-no-sc", "--no-static-condensation",
                                       &static_cond);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_int_option(argc, argv, &it, "-f",
                                      "--flux-averaging", &flux_averaging);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&it, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[it]);
      return 1;
   }

   // 2. Read the beam mesh and verify the expected material and boundary tags.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   dim = CMFEM_Mesh_Dimension(mesh);
   if (CMFEM_Mesh_SpaceDimension(mesh) != dim)
   {
      fprintf(stderr, "invalid mesh\n");
      CMFEM_Mesh_Delete(mesh);
      return 2;
   }
   if (CMFEM_Mesh_AttributesMax(mesh) < 2 ||
       CMFEM_Mesh_BoundaryAttributesMax(mesh) < 2)
   {
      fprintf(stderr,
              "\nInput mesh should have at least two materials and two boundary attributes!\n\n");
      CMFEM_Mesh_Delete(mesh);
      return 3;
   }

   // 3. Convert NURBS meshes to a refined curved mesh that supports AMR.
   if (CMFEM_Mesh_HasNURBSext(mesh))
   {
      CMFEM_Mesh_UniformRefinement(mesh);
      CMFEM_Mesh_UniformRefinement(mesh);
      CMFEM_Mesh_SetCurvature(mesh, 2);
   }

   // 4. Define the vector H1 space for the elasticity displacement field.
   CMFEM_H1_FECollection *fec = CMFEM_H1_FECollection_NewOrderDim(order, dim);
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshH1VDim(
                                          mesh, fec, dim);
   tdim = dim * (dim + 1) / 2;

   // 5. Build the elasticity load vector with the pull-down force on boundary
   //    attribute 2.
   CMFEM_VectorArrayCoefficient *force = CMFEM_VectorArrayCoefficient_New(dim);
   {
      int d;
      for (d = 0; d < dim - 1; d++)
      {
         CMFEM_VectorArrayCoefficient_SetConstantCoefficient(
            force, d, CMFEM_ConstantCoefficient_New(0.0));
      }
   }
   {
      CMFEM_Vector *pull_force = CMFEM_Vector_NewSize(CMFEM_Mesh_BoundaryAttributesMax(
                                                         mesh));
      CMFEM_Vector_Assign(pull_force, 0.0);
      CMFEM_Vector_Set(pull_force, 1, -1.0e-2);
      CMFEM_VectorArrayCoefficient_SetPWConstCoefficient(
         force, dim - 1, CMFEM_PWConstCoefficient_New(pull_force));
      CMFEM_Vector_Delete(pull_force);
   }

   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_LinearForm_AddBoundaryIntegrator_VectorBoundaryLFIntegrator(b, force);

   // 6. Set up the linear elasticity bilinear form and the zero displacement
   //    boundary coefficient used on the clamped side.
   CMFEM_Vector *lambda_values = CMFEM_Vector_NewSize(CMFEM_Mesh_AttributesMax(
                                                         mesh));
   CMFEM_Vector *mu_values = CMFEM_Vector_NewSize(CMFEM_Mesh_AttributesMax(mesh));
   CMFEM_Vector_Assign(lambda_values, 1.0);
   CMFEM_Vector_Assign(mu_values, 1.0);
   CMFEM_Vector_Set(lambda_values, 0, 50.0);
   CMFEM_Vector_Set(mu_values, 0, 50.0);
   CMFEM_PWConstCoefficient *lambda = CMFEM_PWConstCoefficient_New(lambda_values);
   CMFEM_PWConstCoefficient *mu = CMFEM_PWConstCoefficient_New(mu_values);
   CMFEM_Vector_Delete(lambda_values);
   CMFEM_Vector_Delete(mu_values);

   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   CMFEM_BilinearForm_AddDomainIntegrator_Elasticity(a, lambda, mu);
   if (static_cond)
   {
      CMFEM_BilinearForm_EnableStaticCondensation(a);
   }

   CMFEM_Vector *zero_vec = CMFEM_Vector_NewSize(dim);
   CMFEM_Vector_Assign(zero_vec, 0.0);
   CMFEM_VectorConstantCoefficient *zero_vec_coeff =
      CMFEM_VectorConstantCoefficient_New(zero_vec);
   CMFEM_Vector_Delete(zero_vec);

   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_Assign(x, 0.0);

   // 7. Mark boundary attribute 1 as essential.
   _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
      CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
   CMFEM_ArrayInt_Assign(&ess_bdr, 0);
   CMFEM_ArrayInt_Set(&ess_bdr, 0, 1);

   // 8. Build the ZZ estimator on a vector-valued stress recovery space.
   CMFEM_FiniteElementSpace *flux_fespace = CMFEM_FiniteElementSpace_NewMeshH1VDim(
                                               mesh, fec, tdim);
   CMFEM_ElasticityIntegrator *estimator_integ =
      CMFEM_ElasticityIntegrator_NewPWConstCoefficients(lambda, mu);
   CMFEM_ZienkiewiczZhuEstimator *estimator =
      CMFEM_ZienkiewiczZhuEstimator_NewElasticityIntegratorGridFunctionFESpace(
         estimator_integ, x, flux_fespace);
   CMFEM_ZienkiewiczZhuEstimator_SetFluxAveraging(estimator, flux_averaging);

   // 9. Configure the threshold refiner that drives the AMR loop.
   CMFEM_ThresholdRefiner *refiner =
      CMFEM_ThresholdRefiner_NewZienkiewiczZhuEstimator(estimator);
   CMFEM_ThresholdRefiner_SetTotalErrorFraction(refiner, 0.7);

   // 10. Solve on the current mesh, visualize, and adapt until the stopping
   //     criterion or the DOF cap is reached.
   {
      const int max_dofs = 50000;
      const int max_amr_itr = 20;
      for (it = 0; it <= max_amr_itr; it++)
      {
         const int cdofs = CMFEM_FiniteElementSpace_GetTrueVSize(fespace);
         _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list =
            CMFEM_ArrayInt_Construct();
         _Alignas(max_align_t) CMFEM_SparseMatrix A =
            CMFEM_SparseMatrix_Construct();
         _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
         _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();

         printf("\nAMR iteration %d\n", it);
         printf("Number of unknowns: %d\n", cdofs);

         // 11. Assemble the stiffness matrix and load vector on the current
         //     mesh.
         CMFEM_BilinearForm_Assemble(a);
         CMFEM_LinearForm_Assemble(b);

         // 12. Enforce the zero displacement on the essential boundary and
         //     build the constrained linear system.
         CMFEM_GridFunction_ProjectBdrCoefficient_VectorConstantCoefficient(
            x, zero_vec_coeff, &ess_bdr);
         CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr,
                                                       &ess_tdof_list);
         CMFEM_BilinearForm_FormLinearSystemSparseMatrixCopyInterior(
            a, &ess_tdof_list, x, b, &A, &X, &B, 1);

         // 13. Solve the symmetric linear system with PCG and symmetric
         //     Gauss-Seidel preconditioning.
         {
            CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewSparseMatrix(&A);
            CMFEM_PCG_SparseMatrixGSSmoother(&A, M, &B, &X, 3, 2000, 1e-12, 0.0);
            CMFEM_GSSmoother_Delete(M);
         }

         // 14. Recover the displacement field as a finite element grid
         //     function.
         CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);

         // 15. Send the current solution to GLVis.
         if (visualization)
         {
            CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
         }

         CMFEM_SparseMatrix_Destroy(&A);
         CMFEM_Vector_Destroy(&B);
         CMFEM_Vector_Destroy(&X);
         CMFEM_ArrayInt_Destroy(&ess_tdof_list);

         if (cdofs > max_dofs)
         {
            printf("Reached the maximum number of dofs. Stop.\n");
            break;
         }

         // 16. Refine the mesh based on the current error indicators.
         CMFEM_ThresholdRefiner_Apply(refiner, mesh);
         if (CMFEM_ThresholdRefiner_Stop(refiner))
         {
            printf("Stopping criterion satisfied. Stop.\n");
            break;
         }

         // 17. Update the finite element space, the current solution, and the
         //     assembled forms after the mesh changes.
         CMFEM_FiniteElementSpace_Update(fespace);
         CMFEM_GridFunction_Update(x);
         CMFEM_FiniteElementSpace_UpdatesFinished(fespace);
         CMFEM_BilinearForm_Update(a);
         CMFEM_LinearForm_Update(b);
      }
   }

   // 18. Save the reference mesh, a deformed copy of the mesh, and the final
   //     displacement field.
   CMFEM_Mesh_Print(mesh, "ex21_reference.mesh", 16);
   {
      CMFEM_Mesh *deformed = make_deformed_mesh_copy(mesh, fespace, x);
      CMFEM_Mesh_Print(deformed, "ex21_deformed.mesh", 16);
      CMFEM_Mesh_Delete(deformed);
   }
   CMFEM_GridFunction_Save(x, "ex21_displacement.sol", 16);

   // 19. Free the allocated MFEM objects.
   CMFEM_ThresholdRefiner_Delete(refiner);
   CMFEM_ZienkiewiczZhuEstimator_Delete(estimator);
   CMFEM_ElasticityIntegrator_Delete(estimator_integ);
   CMFEM_ArrayInt_Destroy(&ess_bdr);
   CMFEM_GridFunction_Delete(x);
   CMFEM_VectorConstantCoefficient_Delete(zero_vec_coeff);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_LinearForm_Delete(b);
   CMFEM_VectorArrayCoefficient_Delete(force);
   CMFEM_PWConstCoefficient_Delete(lambda);
   CMFEM_PWConstCoefficient_Delete(mu);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_H1_FECollection_Delete(fec);
   CMFEM_Mesh_Delete(mesh);
   return 0;
}
