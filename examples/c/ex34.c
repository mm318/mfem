//                              MFEM C Example 34
//
// Description: This example solves a magnetostatic curl-curl problem where
//              the current density is first computed on a conductor submesh and
//              then transferred back to the full mesh.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("fichera-mixed.mesh");
   int ref_levels = 1;
   int order = 1;
   double delta_const = 1.0e-6;
   int mixed = 1;
   int static_cond = 0;
   int pa = 0;
   const char *device_config = "cpu";
   int visualization = 1;
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
                                         "-mc", "--magnetic-cond",
                                         &delta_const);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-sc", "--static-condensation",
                                       "-no-sc", "--no-static-condensation",
                                       &static_cond);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-mixed", "--mixed-mesh",
                                       "-hex", "--hex-mesh",
                                       &mixed);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-pa", "--partial-assembly",
                                       "-no-pa", "--no-partial-assembly",
                                       &pa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device",
                                         &device_config);
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

   if (!mixed || pa)
   {
      mesh_file = CMFEM_ExamplesDataPath("fichera.mesh");
   }

   // 2. Configure the device and read the mesh.
   {
      CMFEM_Device *device = CMFEM_Device_New(device_config);
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      const int dim = CMFEM_Mesh_Dimension(mesh);
      CMFEM_NdFeCollection *fec_nd = NULL;
      CMFEM_RtFeCollection *fec_rt = NULL;
      CMFEM_FiniteElementSpace *fespace_nd = NULL;
      CMFEM_FiniteElementSpace *fespace_rt = NULL;
      CMFEM_GridFunction *j_full = NULL;
      CMFEM_GridFunction *x = NULL;
      CMFEM_GridFunction *dx = NULL;
      CMFEM_LinearForm *b = NULL;
      CMFEM_BilinearForm *a = NULL;
      CMFEM_DiscreteLinearOperator *curl = NULL;
      CMFEM_ConstantCoefficient *muinv = NULL;
      CMFEM_ConstantCoefficient *delta = NULL;
      _Alignas(max_align_t) CMFEM_ArrayInt submesh_elems = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt sym_plane_attr =
         CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt phi0_attr = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt phi1_attr = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt jn_zero_attr = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
      _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector x_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector dx_vec = CMFEM_Vector_Construct();

      CMFEM_Device_Print(device);

      if (!mixed || pa)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
         if (ref_levels > 0)
         {
            ref_levels--;
         }
      }
      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      if (strcmp(mesh_file, CMFEM_ExamplesDataPath("fichera-mixed.mesh")) == 0)
      {
         CMFEM_ArrayInt_Destroy(&submesh_elems);
         submesh_elems = CMFEM_ArrayInt_ConstructSize(5);
         CMFEM_ArrayInt_Set(&submesh_elems, 0, 0);
         CMFEM_ArrayInt_Set(&submesh_elems, 1, 2);
         CMFEM_ArrayInt_Set(&submesh_elems, 2, 3);
         CMFEM_ArrayInt_Set(&submesh_elems, 3, 4);
         CMFEM_ArrayInt_Set(&submesh_elems, 4, 9);
      }
      else
      {
         CMFEM_ArrayInt_Destroy(&submesh_elems);
         submesh_elems = CMFEM_ArrayInt_ConstructSize(7);
         CMFEM_ArrayInt_Set(&submesh_elems, 0, 10);
         CMFEM_ArrayInt_Set(&submesh_elems, 1, 14);
         CMFEM_ArrayInt_Set(&submesh_elems, 2, 34);
         CMFEM_ArrayInt_Set(&submesh_elems, 3, 36);
         CMFEM_ArrayInt_Set(&submesh_elems, 4, 37);
         CMFEM_ArrayInt_Set(&submesh_elems, 5, 38);
         CMFEM_ArrayInt_Set(&submesh_elems, 6, 39);
      }

      CMFEM_ArrayInt_Destroy(&sym_plane_attr);
      sym_plane_attr = CMFEM_ArrayInt_ConstructSize(8);
      CMFEM_ArrayInt_Set(&sym_plane_attr, 0, 9);
      CMFEM_ArrayInt_Set(&sym_plane_attr, 1, 10);
      CMFEM_ArrayInt_Set(&sym_plane_attr, 2, 11);
      CMFEM_ArrayInt_Set(&sym_plane_attr, 3, 12);
      CMFEM_ArrayInt_Set(&sym_plane_attr, 4, 13);
      CMFEM_ArrayInt_Set(&sym_plane_attr, 5, 14);
      CMFEM_ArrayInt_Set(&sym_plane_attr, 6, 15);
      CMFEM_ArrayInt_Set(&sym_plane_attr, 7, 16);

      CMFEM_ArrayInt_Destroy(&phi0_attr);
      phi0_attr = CMFEM_ArrayInt_ConstructSize(1);
      CMFEM_ArrayInt_Set(&phi0_attr, 0, 2);

      CMFEM_ArrayInt_Destroy(&phi1_attr);
      phi1_attr = CMFEM_ArrayInt_ConstructSize(1);
      CMFEM_ArrayInt_Set(&phi1_attr, 0, 23);

      CMFEM_ArrayInt_Destroy(&jn_zero_attr);
      jn_zero_attr = CMFEM_ArrayInt_ConstructSize(9);
      CMFEM_ArrayInt_Set(&jn_zero_attr, 0, 25);
      for (i = 0; i < CMFEM_ArrayInt_Size(&sym_plane_attr); i++)
      {
         CMFEM_ArrayInt_Set(&jn_zero_attr, i + 1,
                            CMFEM_ArrayInt_Get(&sym_plane_attr, i));
      }

      // 3. Compute the conductor current density on the submesh and transfer
      //    it to the full mesh RT space.
      fec_nd = CMFEM_NdFeCollection_NewOrderDim(order, dim);
      fec_rt = CMFEM_RtFeCollection_NewOrderDim(order - 1, dim);
      fespace_nd = CMFEM_FiniteElementSpace_NewMeshNd(mesh, fec_nd);
      fespace_rt = CMFEM_FiniteElementSpace_NewMeshRt(mesh, fec_rt);
      j_full = CMFEM_GridFunction_New(fespace_rt);
      CMFEM_GridFunction_Assign(j_full, 0.0);
      CMFEM_Ex34_PrepareCurrentDensityRt(mesh,
                                         order,
                                         pa,
                                         visualization,
                                         &submesh_elems,
                                         &phi0_attr,
                                         &phi1_attr,
                                         &jn_zero_attr,
                                         j_full);

      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, j_full, "localhost", 19916);
      }

      // 4. Set up the full H(curl) solve with symmetry-plane essential dofs.
      if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
      {
         CMFEM_ArrayInt_Destroy(&ess_bdr);
         ess_bdr = CMFEM_ArrayInt_ConstructSize(
                      CMFEM_Mesh_BoundaryAttributesMax(mesh));
         CMFEM_ArrayInt_Assign(&ess_bdr, 1);
         for (i = 0; i < CMFEM_ArrayInt_Size(&sym_plane_attr); i++)
         {
            CMFEM_ArrayInt_Set(&ess_bdr,
                               CMFEM_ArrayInt_Get(&sym_plane_attr, i) - 1,
                               0);
         }
         CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace_nd,
                                                       &ess_bdr,
                                                       &ess_tdof_list);
      }

      b = CMFEM_LinearForm_New(fespace_nd);
      CMFEM_Ex34_AssembleCurrentDensityLf(fespace_nd, j_full, b);

      x = CMFEM_GridFunction_New(fespace_nd);
      CMFEM_GridFunction_Assign(x, 0.0);

      muinv = CMFEM_ConstantCoefficient_New(1.0);
      delta = CMFEM_ConstantCoefficient_New(delta_const);
      a = CMFEM_BilinearForm_New(fespace_nd);
      if (pa)
      {
         CMFEM_BilinearForm_SetAssemblyLevelPartial(a);
      }
      CMFEM_BilinearForm_AddDomainIntegratorCci(a, muinv);
      CMFEM_BilinearForm_AddDomainIntegratorVmi(a, delta);
      if (static_cond)
      {
         CMFEM_BilinearForm_EnableStaticCondensation(a);
      }
      CMFEM_BilinearForm_Assemble(a);
      CMFEM_BilinearForm_FormLinearSystemOp(a, &ess_tdof_list, x, b, &A, &X, &B);

      if (pa)
      {
         CMFEM_OperatorJacobiSmoother *M =
            CMFEM_OperatorJacobiSmoother_NewBf(a, &ess_tdof_list);
         printf("\nSolving for magnetic vector potential using CG with a Jacobi preconditioner\n");
         CMFEM_PCGOpOjs(&A, M, &B, &X, 1, 1000, 1.0e-12, 0.0);
         CMFEM_OperatorJacobiSmoother_Delete(M);
      }
      else
      {
         CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewOp(&A);
         printf("\nSolving for magnetic vector potential using CG with a Gauss-Seidel preconditioner\n");
         CMFEM_PCGOpGs(&A, M, &B, &X, 1, 500, 1.0e-12, 0.0);
         CMFEM_GSSmoother_Delete(M);
      }

      CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
      CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
      CMFEM_GridFunction_Save(x, "sol.gf", 8);

      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
      }

      // 5. Compute the magnetic flux as the curl of the vector potential.
      curl = CMFEM_DiscreteLinearOperator_New(fespace_nd, fespace_rt);
      CMFEM_DiscreteLinearOperator_AddDomainInterpolatorCi(curl);
      CMFEM_DiscreteLinearOperator_Assemble(curl);
      dx = CMFEM_GridFunction_New(fespace_rt);
      CMFEM_GridFunction_CopyToVec(x, &x_vec);
      dx_vec = CMFEM_Vector_ConstructSize(CMFEM_FiniteElementSpace_GetVSize(
                                             fespace_rt));
      CMFEM_DiscreteLinearOperator_Mult(curl, &x_vec, &dx_vec);
      CMFEM_GridFunction_SetFromVec(dx, &dx_vec);
      CMFEM_GridFunction_Save(dx, "dsol.gf", 8);

      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, dx, "localhost", 19916);
      }

      CMFEM_Vector_Destroy(&dx_vec);
      CMFEM_Vector_Destroy(&x_vec);
      CMFEM_Vector_Destroy(&B);
      CMFEM_Vector_Destroy(&X);
      CMFEM_OperatorPtr_Destroy(&A);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_ArrayInt_Destroy(&ess_tdof_list);
      CMFEM_ArrayInt_Destroy(&jn_zero_attr);
      CMFEM_ArrayInt_Destroy(&phi1_attr);
      CMFEM_ArrayInt_Destroy(&phi0_attr);
      CMFEM_ArrayInt_Destroy(&sym_plane_attr);
      CMFEM_ArrayInt_Destroy(&submesh_elems);
      CMFEM_DiscreteLinearOperator_Delete(curl);
      CMFEM_ConstantCoefficient_Delete(delta);
      CMFEM_ConstantCoefficient_Delete(muinv);
      CMFEM_BilinearForm_Delete(a);
      CMFEM_GridFunction_Delete(dx);
      CMFEM_GridFunction_Delete(x);
      CMFEM_LinearForm_Delete(b);
      CMFEM_GridFunction_Delete(j_full);
      CMFEM_FiniteElementSpace_Delete(fespace_rt);
      CMFEM_FiniteElementSpace_Delete(fespace_nd);
      CMFEM_RtFeCollection_Delete(fec_rt);
      CMFEM_NdFeCollection_Delete(fec_nd);
      CMFEM_Mesh_Delete(mesh);
      CMFEM_Device_Delete(device);
   }

   return 0;
}
