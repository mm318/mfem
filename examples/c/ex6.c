//                                MFEM C Example 6
//
// Description: This is a version of Example 1 with a simple adaptive mesh
//              refinement loop. The problem being solved is the Poisson
//              equation -Delta u = 1 with homogeneous Dirichlet boundary
//              conditions.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int pa = 0;
   const char *device_config = "cpu";
   int max_dofs = 50000;
   int lszz = 0;
   int visualization = 1;
   int dim;
   int sdim;
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
                                       "-pa", "--partial-assembly",
                                       "-no-pa", "--no-partial-assembly",
                                       &pa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_string_option(argc, argv, &it, "-d", "--device",
                                         &device_config);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &it, "-md", "--max-dofs",
                                      &max_dofs);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&it, argv,
                                       "-ls", "--ls-zz",
                                       "-no-ls", "--no-ls-zz",
                                       &lszz);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&it, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[it]);
      return 1;
   }

   // 2. Enable hardware devices and programming models.
   {
      CMFEM_Device *device = CMFEM_Device_New(device_config);
      CMFEM_Device_Print(device);

      // 3. Read the mesh from the given mesh file.
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      dim = CMFEM_Mesh_Dimension(mesh);
      sdim = CMFEM_Mesh_SpaceDimension(mesh);

      // 4. Convert NURBS meshes to piecewise-polynomial curved meshes.
      if (CMFEM_Mesh_HasNURBSext(mesh))
      {
         for (it = 0; it < 2; it++)
         {
            CMFEM_Mesh_UniformRefinement(mesh);
         }
         CMFEM_Mesh_SetCurvature(mesh, 2);
      }

      // 5. Define a finite element space on the mesh.
      CMFEM_H1FeCollection *fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshH1(
                                             mesh, fec);

      // 6. Set up the bilinear and linear forms.
      CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
      if (pa)
      {
         CMFEM_BilinearForm_SetAssemblyLevelPartial(a);
         CMFEM_BilinearForm_SetDiagonalPolicyOne(a);
      }
      CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);

      CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
      CMFEM_ConstantCoefficient *zero = CMFEM_ConstantCoefficient_New(0.0);
      CMFEM_BilinearForm_AddDomainIntegratorDiCc(a, one);
      CMFEM_LinearForm_AddDomainIntegratorDliCc(
         b, one);

      // 7. Initialize the solution grid function.
      CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
      CMFEM_GridFunction_Assign(x, 0.0);

      // 8. Mark all boundary attributes as essential.
      if (CMFEM_Mesh_BoundaryAttributesSize(mesh) <= 0)
      {
         fprintf(stderr, "Boundary attributes required in the mesh.\n");
         CMFEM_GridFunction_Delete(x);
         CMFEM_ConstantCoefficient_Delete(zero);
         CMFEM_ConstantCoefficient_Delete(one);
         CMFEM_LinearForm_Delete(b);
         CMFEM_BilinearForm_Delete(a);
         CMFEM_FiniteElementSpace_Delete(fespace);
         CMFEM_H1FeCollection_Delete(fec);
         CMFEM_Mesh_Delete(mesh);
         CMFEM_Device_Delete(device);
         return 1;
      }
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
         CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr, 1);

      // 9. Set up an error estimator.
      CMFEM_DiffusionIntegrator *diff_integ =
         CMFEM_DiffusionIntegrator_NewCc(one);
      CMFEM_ThresholdRefiner *refiner = NULL;
      CMFEM_ZienkiewiczZhuEstimator *zz_estimator = NULL;
      CMFEM_LSZienkiewiczZhuEstimator *lszz_estimator = NULL;

      if (lszz)
      {
         lszz_estimator =
            CMFEM_LSZienkiewiczZhuEstimator_NewDiGf(
               diff_integ, x);
         if (dim == 3 &&
             CMFEM_Mesh_GetElementType(mesh, 0) != CMFEM_ELEMENT_HEXAHEDRON)
         {
            CMFEM_LSZienkiewiczZhuEstimator_SetTichonovRegularization(
               lszz_estimator, 1.0e-8);
         }
         refiner = CMFEM_ThresholdRefiner_NewLzz(
                      lszz_estimator);
      }
      else
      {
         CMFEM_FiniteElementSpace *flux_fes =
            CMFEM_FiniteElementSpace_NewMeshH1VDim(mesh, fec, sdim);
         zz_estimator =
            CMFEM_ZienkiewiczZhuEstimator_NewDiGfFes(
               diff_integ, x, flux_fes);
         CMFEM_ZienkiewiczZhuEstimator_SetAnisotropic(zz_estimator, 1);
         refiner = CMFEM_ThresholdRefiner_NewZze(
                      zz_estimator);
      }
      CMFEM_ThresholdRefiner_SetTotalErrorFraction(refiner, 0.7);

      // 10. Main AMR loop.
      for (it = 0; ; it++)
      {
         const int cdofs = CMFEM_FiniteElementSpace_GetTrueVSize(fespace);
         printf("\nAMR iteration %d\n", it);
         printf("Number of unknowns: %d\n", cdofs);

         // 11. Assemble the right-hand side and set Dirichlet values.
         CMFEM_LinearForm_Assemble(b);
         CMFEM_GridFunction_ProjectBdrCoefficientCc(x, zero,
                                                    &ess_bdr);
         _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
         CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr,
                                                       &ess_tdof_list);

         // 12. Assemble the stiffness matrix and form the linear system.
         CMFEM_BilinearForm_Assemble(a);
         _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
         _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
         _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
         CMFEM_BilinearForm_FormLinearSystemOpCopyInterior(a,
                                                           &ess_tdof_list,
                                                           x,
                                                           b,
                                                           &A,
                                                           &X,
                                                           &B,
                                                           1);

         // 13. Solve the linear system.
         if (!pa)
         {
            CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewOp(&A);
            CMFEM_PCGOpGs(&A, M, &B, &X, 3, 200, 1.0e-12, 0.0);
            CMFEM_GSSmoother_Delete(M);
         }
         else
         {
            CMFEM_OperatorJacobiSmoother *M =
               CMFEM_OperatorJacobiSmoother_NewBf(a, &ess_tdof_list);
            CMFEM_PCGOpOjs(&A, M, &B, &X, 3, 2000, 1.0e-12,
                           0.0);
            CMFEM_OperatorJacobiSmoother_Delete(M);
         }

         // 14. Recover the finite element solution and visualize it.
         CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
         if (visualization)
         {
            CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
         }

         CMFEM_OperatorPtr_Destroy(&A);
         CMFEM_Vector_Destroy(&B);
         CMFEM_Vector_Destroy(&X);
         CMFEM_ArrayInt_Destroy(&ess_tdof_list);

         if (cdofs > max_dofs)
         {
            printf("Reached the maximum number of dofs. Stop.\n");
            break;
         }

         // 15. Refine the mesh and update the dependent objects.
         if (!CMFEM_ThresholdRefiner_Apply(refiner, mesh))
         {
            if (CMFEM_ThresholdRefiner_Stop(refiner))
            {
               printf("Stopping criterion satisfied. Stop.\n");
            }
            break;
         }

         CMFEM_FiniteElementSpace_Update(fespace);
         CMFEM_GridFunction_Update(x);
         CMFEM_BilinearForm_Update(a);
         CMFEM_LinearForm_Update(b);
      }

      CMFEM_ThresholdRefiner_Delete(refiner);
      if (zz_estimator)
      {
         CMFEM_ZienkiewiczZhuEstimator_Delete(zz_estimator);
      }
      if (lszz_estimator)
      {
         CMFEM_LSZienkiewiczZhuEstimator_Delete(lszz_estimator);
      }
      CMFEM_DiffusionIntegrator_Delete(diff_integ);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_GridFunction_Delete(x);
      CMFEM_ConstantCoefficient_Delete(zero);
      CMFEM_ConstantCoefficient_Delete(one);
      CMFEM_LinearForm_Delete(b);
      CMFEM_BilinearForm_Delete(a);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_H1FeCollection_Delete(fec);
      CMFEM_Mesh_Delete(mesh);
      CMFEM_Device_Delete(device);
   }

   return 0;
}
