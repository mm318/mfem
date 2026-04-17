//                              MFEM C Example 39
//
// Description: This example code demonstrates named attribute sets in MFEM
//              and uses them in a simple Poisson solve with inhomogeneous
//              diffusion coefficients.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("compass.msh");
   int order = 1;
   const char *source_name = "Rose Even";
   const char *ess_name = "Boundary";
   int visualization = 1;
   int dim;
   int ref_levels;
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

      parsed = cmfem_parse_string_option(argc, argv, &i, "-src",
                                         "--source-attr-name", &source_name);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-ess",
                                         "--ess-attr-name", &ess_name);
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

   // 2. Read the mesh from the given mesh file.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   dim = CMFEM_Mesh_Dimension(mesh);

   // 3. Refine the mesh to increase the resolution.
   ref_levels = cmfem_uniform_refinement_levels(50000.0, CMFEM_Mesh_GetNE(mesh),
                                                dim);
   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }

   // 4. Display the named attribute sets and build the derived sets used by
   //    this example.
   CMFEM_Mesh_PrintAttributeSetNames(mesh, 0, "Element Attribute Set Names");
   CMFEM_Mesh_PrintAttributeSetNames(mesh, 1, "Boundary Attribute Set Names");

   CMFEM_Mesh_SetAttributeSetFromName(mesh, 0, "North", "N Even");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "North", "N Odd");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 0, "South", "S Even");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "South", "S Odd");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 0, "East", "E Even");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "East", "E Odd");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 0, "West", "W Even");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "West", "W Odd");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 0, "Rose Even", "N Even");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "Rose Even", "S Even");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "Rose Even", "E Even");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "Rose Even", "W Even");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 0, "Rose Odd", "N Odd");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "Rose Odd", "S Odd");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "Rose Odd", "E Odd");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "Rose Odd", "W Odd");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 0, "Rose", "Rose Even");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 0, "Rose", "Rose Odd");

   CMFEM_Mesh_SetAttributeSetFromName(mesh, 1, "Northern Boundary", "NNE");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 1, "Northern Boundary", "NNW");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 1, "Southern Boundary", "SSE");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 1, "Southern Boundary", "SSW");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 1, "Eastern Boundary", "ENE");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 1, "Eastern Boundary", "ESE");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 1, "Western Boundary", "WNW");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 1, "Western Boundary", "WSW");
   CMFEM_Mesh_SetAttributeSetFromName(mesh, 1, "Boundary", "Northern Boundary");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 1, "Boundary", "Southern Boundary");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 1, "Boundary", "Eastern Boundary");
   CMFEM_Mesh_AddToAttributeSetFromName(mesh, 1, "Boundary", "Western Boundary");

   // 5. Define a finite element space on the mesh.
   CMFEM_H1FeCollection *fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh,
                                                                          fec);
   printf("Number of finite element unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   // 6. Determine the list of true essential boundary dofs.
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
   if (CMFEM_Mesh_AttributeSetExists(mesh, 1, ess_name))
   {
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr_marker =
         CMFEM_ArrayInt_Construct();
      CMFEM_Mesh_GetAttributeSetMarker(mesh, 1, ess_name, &ess_bdr_marker);
      CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr_marker,
                                                    &ess_tdof_list);
      CMFEM_ArrayInt_Destroy(&ess_bdr_marker);
   }

   // 7. Set up the linear form b(.) corresponding to a source indicator
   //    function over the named attribute set `source_name`.
   _Alignas(max_align_t) CMFEM_ArrayInt source_marker = CMFEM_ArrayInt_Construct();
   CMFEM_Mesh_GetAttributeSetMarker(mesh, 0, source_name, &source_marker);
   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_LinearForm_AddDomainIntegratorDliCcMarker(
      b, one, &source_marker);
   CMFEM_LinearForm_Assemble(b);

   // 8. Define the solution vector x corresponding to fespace and initialize it
   //    with zero.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_Assign(x, 0.0);

   // 9. Set up the bilinear form with default, base, and rose-region diffusion
   //    coefficients.
   _Alignas(max_align_t) CMFEM_ArrayInt base_marker = CMFEM_ArrayInt_Construct();
   _Alignas(max_align_t) CMFEM_ArrayInt rose_marker = CMFEM_ArrayInt_Construct();
   CMFEM_Mesh_GetAttributeSetMarker(mesh, 0, "Base", &base_marker);
   CMFEM_Mesh_GetAttributeSetMarker(mesh, 0, "Rose Even", &rose_marker);

   CMFEM_ConstantCoefficient *default_coef =
      CMFEM_ConstantCoefficient_New(1.0e-6);
   CMFEM_ConstantCoefficient *base_coef = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_ConstantCoefficient *rose_coef = CMFEM_ConstantCoefficient_New(2.0);

   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   CMFEM_BilinearForm_AddDomainIntegratorDiCc(a, default_coef);
   CMFEM_BilinearForm_AddDomainIntegratorDiCcMarker(
      a, base_coef, &base_marker);
   CMFEM_BilinearForm_AddDomainIntegratorDiCcMarker(
      a, rose_coef, &rose_marker);
   CMFEM_BilinearForm_Assemble(a);

   // 10. Assemble and solve the linear system.
   _Alignas(max_align_t) CMFEM_SparseMatrix A = CMFEM_SparseMatrix_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
   CMFEM_BilinearForm_FormLinearSystemSm(a, &ess_tdof_list, x, b, &A, &X,
                                         &B);
   printf("Size of linear system: %d\n", CMFEM_SparseMatrix_Height(&A));

   {
      CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewSm(&A);
      CMFEM_PCGSmGs(&A, M, &B, &X, 1, 800, 1e-12, 0.0);
      CMFEM_GSSmoother_Delete(M);
   }

   // 11. Recover the solution as a finite element grid function.
   CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);

   // 12. Save the refined mesh and the solution.
   CMFEM_Mesh_Save(mesh, "refined.mesh", 8);
   CMFEM_GridFunction_Save(x, "sol.gf", 8);

   // 13. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
   }

   CMFEM_SparseMatrix_Destroy(&A);
   CMFEM_Vector_Destroy(&B);
   CMFEM_Vector_Destroy(&X);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_ConstantCoefficient_Delete(rose_coef);
   CMFEM_ConstantCoefficient_Delete(base_coef);
   CMFEM_ConstantCoefficient_Delete(default_coef);
   CMFEM_GridFunction_Delete(x);
   CMFEM_LinearForm_Delete(b);
   CMFEM_ConstantCoefficient_Delete(one);
   CMFEM_ArrayInt_Destroy(&rose_marker);
   CMFEM_ArrayInt_Destroy(&base_marker);
   CMFEM_ArrayInt_Destroy(&source_marker);
   CMFEM_ArrayInt_Destroy(&ess_tdof_list);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_H1FeCollection_Delete(fec);
   CMFEM_Mesh_Delete(mesh);
   return 0;
}
