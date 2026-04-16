const std = @import("std");

const Precision = enum {
    double,
    single,
};

const general_sources = &.{
    "core/general/array.cpp",
    "core/general/binaryio.cpp",
    "core/general/cuda.cpp",
    "core/general/device.cpp",
    "core/general/error.cpp",
    "core/general/gecko.cpp",
    "core/general/globals.cpp",
    "core/general/hash.cpp",
    "core/general/hash_util.cpp",
    "core/general/isockstream.cpp",
    "core/general/mem_manager.cpp",
    "core/general/occa.cpp",
    "core/general/optparser.cpp",
    "core/general/osockstream.cpp",
    "core/general/sets.cpp",
    "core/general/socketstream.cpp",
    "core/general/stable3d.cpp",
    "core/general/table.cpp",
    "core/general/tic_toc.cpp",
    "core/general/tinyxml2.cpp",
    "core/general/version.cpp",
    "core/general/hip.cpp",
};

const linalg_sources = &.{
    "core/linalg/auxiliary.cpp",
    "core/linalg/batched/batched.cpp",
    "core/linalg/batched/gpu_blas.cpp",
    "core/linalg/batched/magma.cpp",
    "core/linalg/batched/native.cpp",
    "core/linalg/batched/solver.cpp",
    "core/linalg/blockmatrix.cpp",
    "core/linalg/blockoperator.cpp",
    "core/linalg/blockvector.cpp",
    "core/linalg/complex_densemat.cpp",
    "core/linalg/complex_operator.cpp",
    "core/linalg/constraints.cpp",
    "core/linalg/densemat.cpp",
    "core/linalg/filteredsolver.cpp",
    "core/linalg/handle.cpp",
    "core/linalg/matrix.cpp",
    "core/linalg/mma.cpp",
    "core/linalg/ode.cpp",
    "core/linalg/operator.cpp",
    "core/linalg/ordering.cpp",
    "core/linalg/particlevector.cpp",
    "core/linalg/solvers.cpp",
    "core/linalg/sparsemat.cpp",
    "core/linalg/sparsesmoothers.cpp",
    "core/linalg/symmat.cpp",
    "core/linalg/vector.cpp",
};

const mesh_sources = &.{
    "core/mesh/attribute_sets.cpp",
    "core/mesh/element.cpp",
    "core/mesh/exodus_writer.cpp",
    "core/mesh/face_nbr_geom.cpp",
    "core/mesh/gmsh.cpp",
    "core/mesh/hexahedron.cpp",
    "core/mesh/mesh.cpp",
    "core/mesh/mesh_operators.cpp",
    "core/mesh/mesh_readers.cpp",
    "core/mesh/ncmesh.cpp",
    "core/mesh/ncnurbs.cpp",
    "core/mesh/nurbs.cpp",
    "core/mesh/point.cpp",
    "core/mesh/pyramid.cpp",
    "core/mesh/quadrilateral.cpp",
    "core/mesh/segment.cpp",
    "core/mesh/spacing.cpp",
    "core/mesh/tetrahedron.cpp",
    "core/mesh/triangle.cpp",
    "core/mesh/vertex.cpp",
    "core/mesh/vtk.cpp",
    "core/mesh/vtkhdf.cpp",
    "core/mesh/wedge.cpp",
    "core/mesh/submesh/submesh.cpp",
    "core/mesh/submesh/ncsubmesh.cpp",
    "core/mesh/submesh/submesh_utils.cpp",
    "core/mesh/submesh/transfermap.cpp",
};

const fem_sources = &.{
    "core/fem/bilinearform.cpp",
    "core/fem/bilinearform_ext.cpp",
    "core/fem/bilininteg.cpp",
    "core/fem/integ/bilininteg_br2.cpp",
    "core/fem/integ/bilininteg_convection_mf.cpp",
    "core/fem/integ/bilininteg_convection_pa.cpp",
    "core/fem/integ/bilininteg_convection_ea.cpp",
    "core/fem/integ/bilininteg_curlcurl_pa.cpp",
    "core/fem/integ/bilininteg_dgdiffusion_pa.cpp",
    "core/fem/integ/bilininteg_dgtrace_pa.cpp",
    "core/fem/integ/bilininteg_dgtrace_ea.cpp",
    "core/fem/integ/bilininteg_diffusion_mf.cpp",
    "core/fem/integ/bilininteg_diffusion_pa.cpp",
    "core/fem/integ/bilininteg_diffusion_ea.cpp",
    "core/fem/integ/bilininteg_diffusion_patch.cpp",
    "core/fem/integ/bilininteg_divdiv_pa.cpp",
    "core/fem/integ/bilininteg_elasticity_ea.cpp",
    "core/fem/integ/bilininteg_elasticity_pa.cpp",
    "core/fem/integ/bilininteg_gradient_pa.cpp",
    "core/fem/integ/bilininteg_interp_pa.cpp",
    "core/fem/integ/bilininteg_mass_mf.cpp",
    "core/fem/integ/bilininteg_mass_pa.cpp",
    "core/fem/integ/bilininteg_mass_ea.cpp",
    "core/fem/integ/bilininteg_mixedcurl_pa.cpp",
    "core/fem/integ/bilininteg_mixedvecgrad_pa.cpp",
    "core/fem/integ/bilininteg_trace_jump_ea.cpp",
    "core/fem/integ/bilininteg_transpose_ea.cpp",
    "core/fem/integ/bilininteg_vecdiffusion_mf.cpp",
    "core/fem/integ/bilininteg_vecdiffusion_pa.cpp",
    "core/fem/integ/bilininteg_vecdiv_pa.cpp",
    "core/fem/integ/bilininteg_vecmass_mf.cpp",
    "core/fem/integ/bilininteg_vecmass_pa.cpp",
    "core/fem/integ/bilininteg_vectorfediv_pa.cpp",
    "core/fem/integ/bilininteg_vectorfemass_pa.cpp",
    "core/fem/integ/bilininteg_diffusion_kernels.cpp",
    "core/fem/integ/bilininteg_elasticity_kernels.cpp",
    "core/fem/integ/bilininteg_hcurl_kernels.cpp",
    "core/fem/integ/bilininteg_hdiv_ea.cpp",
    "core/fem/integ/bilininteg_hdiv_kernels.cpp",
    "core/fem/integ/bilininteg_hcurlhdiv_kernels.cpp",
    "core/fem/integ/bilininteg_mass_kernels.cpp",
    "core/fem/integ/lininteg_boundary.cpp",
    "core/fem/integ/lininteg_boundary_flux.cpp",
    "core/fem/integ/lininteg_domain.cpp",
    "core/fem/integ/lininteg_domain_grad.cpp",
    "core/fem/integ/lininteg_domain_vectorfe.cpp",
    "core/fem/integ/nonlininteg_vecconvection_pa.cpp",
    "core/fem/integ/nonlininteg_vecconvection_mf.cpp",
    "core/fem/coefficient.cpp",
    "core/fem/complex_fem.cpp",
    "core/fem/convergence.cpp",
    "core/fem/datacollection.cpp",
    "core/fem/dgmassinv.cpp",
    "core/fem/doftrans.cpp",
    "core/fem/dfem/doperator.cpp",
    "core/fem/eltrans.cpp",
    "core/fem/batchitrans.cpp",
    "core/fem/estimators.cpp",
    "core/fem/fe.cpp",
    "core/fem/fe/face_map_utils.cpp",
    "core/fem/fe/fe_base.cpp",
    "core/fem/fe/fe_fixed_order.cpp",
    "core/fem/fe/fe_h1.cpp",
    "core/fem/fe/fe_l2.cpp",
    "core/fem/fe/fe_nd.cpp",
    "core/fem/fe/fe_nurbs.cpp",
    "core/fem/fe/fe_pos.cpp",
    "core/fem/fe/fe_pyramid.cpp",
    "core/fem/fe/fe_rt.cpp",
    "core/fem/fe/fe_ser.cpp",
    "core/fem/fe_coll.cpp",
    "core/fem/fespace.cpp",
    "core/fem/derefmat_op.cpp",
    "core/fem/pderefmat_op.cpp",
    "core/fem/geom.cpp",
    "core/fem/gridfunc.cpp",
    "core/fem/hybridization.cpp",
    "core/fem/hybridization_ext.cpp",
    "core/fem/intrules.cpp",
    "core/fem/intrules_cut.cpp",
    "core/fem/ceed/interface/basis.cpp",
    "core/fem/ceed/interface/restriction.cpp",
    "core/fem/ceed/interface/operator.cpp",
    "core/fem/ceed/interface/util.cpp",
    "core/fem/ceed/integrators/convection/convection.cpp",
    "core/fem/ceed/integrators/diffusion/diffusion.cpp",
    "core/fem/ceed/integrators/nlconvection/nlconvection.cpp",
    "core/fem/ceed/integrators/mass/mass.cpp",
    "core/fem/ceed/solvers/algebraic.cpp",
    "core/fem/ceed/solvers/full-assembly.cpp",
    "core/fem/ceed/solvers/solvers-atpmg.cpp",
    "core/fem/kdtree.cpp",
    "core/fem/linearform.cpp",
    "core/fem/linearform_ext.cpp",
    "core/fem/lininteg.cpp",
    "core/fem/lor/lor.cpp",
    "core/fem/lor/lor_ads.cpp",
    "core/fem/lor/lor_ams.cpp",
    "core/fem/lor/lor_batched.cpp",
    "core/fem/multigrid.cpp",
    "core/fem/nonlinearform.cpp",
    "core/fem/nonlinearform_ext.cpp",
    "core/fem/nonlininteg.cpp",
    "core/fem/fespacehierarchy.cpp",
    "core/fem/qfunction.cpp",
    "core/fem/qinterp/det.cpp",
    "core/fem/qinterp/eval_by_nodes.cpp",
    "core/fem/qinterp/eval_by_vdim.cpp",
    "core/fem/qinterp/eval_hdiv.cpp",
    "core/fem/qinterp/grad_by_nodes.cpp",
    "core/fem/qinterp/grad_by_vdim.cpp",
    "core/fem/qspace.cpp",
    "core/fem/quadinterpolator.cpp",
    "core/fem/quadinterpolator_face.cpp",
    "core/fem/restriction.cpp",
    "core/fem/normal_deriv_restriction.cpp",
    "core/fem/staticcond.cpp",
    "core/fem/tmop.cpp",
    "core/fem/tmop/pa.cpp",
    "core/fem/tmop/assemble/diag2_limit.cpp",
    "core/fem/tmop/assemble/diag2.cpp",
    "core/fem/tmop/assemble/grad2_limit.cpp",
    "core/fem/tmop/assemble/grad2.cpp",
    "core/fem/tmop/assemble/diag3_limit.cpp",
    "core/fem/tmop/assemble/diag3.cpp",
    "core/fem/tmop/assemble/grad3_limit.cpp",
    "core/fem/tmop/assemble/grad3.cpp",
    "core/fem/tmop/metrics/001.cpp",
    "core/fem/tmop/metrics/002.cpp",
    "core/fem/tmop/metrics/007.cpp",
    "core/fem/tmop/metrics/056.cpp",
    "core/fem/tmop/metrics/077.cpp",
    "core/fem/tmop/metrics/080.cpp",
    "core/fem/tmop/metrics/094.cpp",
    "core/fem/tmop/metrics/302.cpp",
    "core/fem/tmop/metrics/303.cpp",
    "core/fem/tmop/metrics/315.cpp",
    "core/fem/tmop/metrics/318.cpp",
    "core/fem/tmop/metrics/321.cpp",
    "core/fem/tmop/metrics/332.cpp",
    "core/fem/tmop/metrics/338.cpp",
    "core/fem/tmop/mult/grad2_limit.cpp",
    "core/fem/tmop/mult/grad2.cpp",
    "core/fem/tmop/mult/mult2_limit.cpp",
    "core/fem/tmop/mult/mult2.cpp",
    "core/fem/tmop/mult/grad3_limit.cpp",
    "core/fem/tmop/mult/grad3.cpp",
    "core/fem/tmop/mult/mult3_limit.cpp",
    "core/fem/tmop/mult/mult3.cpp",
    "core/fem/tmop/tools/det2_jpr.cpp",
    "core/fem/tmop/tools/det3_jpr.cpp",
    "core/fem/tmop/tools/discrete.cpp",
    "core/fem/tmop/tools/energy2_limit.cpp",
    "core/fem/tmop/tools/energy2.cpp",
    "core/fem/tmop/tools/energy3_limit.cpp",
    "core/fem/tmop/tools/energy3.cpp",
    "core/fem/tmop/tools/target2.cpp",
    "core/fem/tmop/tools/target3.cpp",
    "core/fem/tmop_tools.cpp",
    "core/fem/tmop_amr.cpp",
    "core/fem/gslib.cpp",
    "core/fem/gslib/findpts_local_2.cpp",
    "core/fem/gslib/findpts_local_3.cpp",
    "core/fem/gslib/interpolate_local_2.cpp",
    "core/fem/gslib/interpolate_local_3.cpp",
    "core/fem/transfer.cpp",
    "core/fem/hyperbolic.cpp",
    "core/fem/integrator.cpp",
    "core/fem/bounds.cpp",
    "core/fem/particleset.cpp",
};

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const shared = b.option(bool, "shared", "Build MFEM as a shared library") orelse false;
    const precision = b.option(Precision, "precision", "Floating-point precision") orelse .double;
    const use_exceptions = b.option(bool, "exceptions", "Enable MFEM exceptions") orelse false;
    const thread_safe = b.option(bool, "thread-safe", "Enable MFEM thread safety") orelse false;
    const use_memalloc = b.option(bool, "memalloc", "Enable MFEM internal MEMALLOC") orelse true;
    const use_simd = b.option(bool, "simd", "Enable MFEM SIMD code paths") orelse false;

    const write_files = b.addWriteFiles();
    _ = write_files.add("config/_config.hpp", mfemConfigHeader(
        b,
        target,
        optimize,
        shared,
        precision,
        use_exceptions,
        thread_safe,
        use_memalloc,
        use_simd,
    ));

    const mfem = b.addLibrary(.{
        .name = "mfem",
        .linkage = if (shared) .dynamic else .static,
        .root_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libc = true,
            .link_libcpp = true,
        }),
    });

    mfem.root_module.addIncludePath(write_files.getDirectory());
    mfem.root_module.addIncludePath(b.path("config"));
    mfem.root_module.addIncludePath(b.path("."));
    mfem.root_module.addCMacro("MFEM_CONFIG_FILE", "\"config/_config.hpp\"");

    if (shared and target.result.os.tag == .windows) {
        mfem.root_module.addCMacro("mfem_EXPORTS", "1");
    }

    addSourceGroup(mfem.root_module, b, general_sources);
    addSourceGroup(mfem.root_module, b, linalg_sources);
    addSourceGroup(mfem.root_module, b, mesh_sources);
    addSourceGroup(mfem.root_module, b, fem_sources);

    b.installArtifact(mfem);
}

fn addSourceGroup(
    module: *std.Build.Module,
    b: *std.Build,
    files: []const []const u8,
) void {
    module.addCSourceFiles(.{
        .root = b.path("."),
        .files = files,
        .flags = &.{
            "-std=c++17",
        },
    });
}

fn mfemConfigHeader(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    shared: bool,
    precision: Precision,
    use_exceptions: bool,
    thread_safe: bool,
    use_memalloc: bool,
    use_simd: bool,
) []const u8 {
    const timer_type: u8 = switch (target.result.os.tag) {
        .windows => 3,
        .macos, .ios, .tvos, .watchos, .visionos => 4,
        else => 2,
    };

    return b.fmt(
        \\#ifndef MFEM_CONFIG_HEADER
        \\#define MFEM_CONFIG_HEADER
        \\
        \\#define MFEM_VERSION 40901
        \\#define MFEM_VERSION_STRING "4.9.1"
        \\#define MFEM_VERSION_MAJOR ((MFEM_VERSION) / 10000)
        \\#define MFEM_VERSION_MINOR (((MFEM_VERSION) / 100) % 100)
        \\#define MFEM_VERSION_PATCH ((MFEM_VERSION) % 100)
        \\#define MFEM_SOURCE_DIR "."
        \\#define MFEM_INSTALL_DIR "zig-out"
        \\#define MFEM_GIT_STRING "(unknown)"
        \\{s}
        \\{s}
        \\{s}
        \\{s}
        \\{s}
        \\{s}
        \\{s}
        \\#define MFEM_TIMER_TYPE {d}
        \\
        \\#endif // MFEM_CONFIG_HEADER
        \\
    ,
        .{
            if (shared) "#define MFEM_SHARED_BUILD" else "",
            switch (precision) {
                .double => "#define MFEM_USE_DOUBLE",
                .single => "#define MFEM_USE_SINGLE",
            },
            if (optimize == .Debug) "#define MFEM_DEBUG" else "",
            if (use_exceptions) "#define MFEM_USE_EXCEPTIONS" else "",
            if (thread_safe) "#define MFEM_THREAD_SAFE" else "",
            if (use_memalloc) "#define MFEM_USE_MEMALLOC" else "",
            if (use_simd) "#define MFEM_USE_SIMD" else "",
            timer_type,
        },
    );
}
