# MFEM Agent Notes

This file supplements `README.md`, `INSTALL`, and `CONTRIBUTING.md`. It is for
repo-specific conventions and subtle architecture/build facts that are easy to
miss when editing `core/`, `c_api/`, and the Zig build.

## Git Hooks

- Recommended git hooks live under `config/githooks/`.
- `config/githooks/pre-commit` runs `astyle` with the repo style settings.
- `config/githooks/pre-push` runs the quick repository checks used by CI, and
`./config/githooks/pre-push --style` is the fastest way to reproduce the
GitHub Actions `repo-check.yml` style failure locally.
- When code style changes are involved, prefer validating with the hook itself
instead of only running `astyle` manually, because the hook is what the repo
and workflow actually enforce.

## C API Wrapper Rules

- The public C entry header is `c_api/cmfem.h`. Keep it as an aggregate header
only;
do not put implementations or shared macros there.
- Put shared public C declarations in `c_api/common.h`.
- Put shared internal C++-only wrapper helpers in `c_api/common.hpp`.
- Mirror MFEM's subsystem layout under `c_api/`.
  Put wrappers for core MFEM types under paths such as
  `c_api/general/...`, `c_api/linalg/...`, `c_api/mesh/...`, and
  `c_api/fem/...`, matching the subsystem where the wrapped C++ type lives.
- Example-specific focused helper bridge headers and implementations belong
  under `examples/c/adapters/exNN/...`, matching the example/helper source
  they wrap.
- Those adapter headers should depend on `c_api/` headers such as
  `c_api/common.h`; `c_api/` must not depend on the adapter layer.
- Each wrapped MFEM type still gets its own header/implementation pair. Do not
  add new wrappers back into a monolithic file.
- The public storage type shape is fixed:
`struct CMFEM_<Name> { char data[N]; }`.
- Use `CMFEM_STORAGE(Name, Size)` in the per-type public header to define the
   storage struct. Do not add hidden pointers, alignment members, or ownership
   flags to the public C struct.
   - Normalize wrapper class/type names to a single `CMFEM_<ClassName>` segment.
      Do not put extra underscores inside the class/type segment. For example, use
         `CMFEM_H1FeCollection`, `CMFEM_DgFeCollection`, and `CMFEM_ArrayInt`, not
         names like `CMFEM_H1_FECollection`.
         - If the wrapped C++ type is a template specialization, bake the substituted
         type into the C name using a compact 1-3 character PascalCase shorthand. For
         example, `mfem::Array<int>` maps to `CMFEM_ArrayInt`.
- In the `.cpp` implementation, `static_assert` the storage size with
   `CMFEM_ASSERT_TYPE(CMFEM_<Name>, mfem::<CppType>)`.
   - Reuse the casting helpers in `c_api/common.hpp`, especially
   `cmfem::As<T>()`. Do not hand-roll `reinterpret_cast` helpers in each file.
   - Reuse the shared reference helpers in `c_api/common.hpp` when they already
   exist, such as `ArrayIntRef`, `VectorRef`, `OperatorPtrRef`, and
   `SparseMatrixFromOperator`.
   - Member-function wrappers are exposed as free C functions named
   `CMFEM_<ClassName>_<MethodName>(...)`.
   - Free MFEM functions stay free C functions and are named
   `CMFEM_<FunctionName>(...)`.
   - The only underscores in public wrapper function names are the separator after
   `CMFEM` and, for member wrappers, the separator between the class segment and
      the method segment. Do not put extra underscores inside the method segment.
      - If a wrapper name needs to encode template substitutions or overload
      disambiguation types, append those type shorthands directly to the class or
      method segment with no extra underscores. Use compact 1-3 character
      PascalCase shorthands such as:
      `Cc` for `ConstantCoefficient`, `Fc` for `FunctionCoefficient`,
            `Gfc` for `GridFunctionCoefficient`, `Vcc` for `VectorConstantCoefficient`,
                  `Vfc` for
                     `VectorFunctionCoefficient`, `Pwc` for `PWConstCoefficient`, `Ai` for
                           `ArrayInt`, `Sm` for `SparseMatrix`, `Op` for
                                 `OperatorPtr`, `Gf` for `GridFunction`, `Fes` for `FiniteElementSpace`,
                                       `Di` for `DiffusionIntegrator`, `Ei` for `ElasticityIntegrator`,
                                             `Bli` for `BoundaryLFIntegrator`, `Bmi` for
                                                   `BoundaryMassIntegrator`,
                                                   `Zze` for `ZienkiewiczZhuEstimator`, `Lzz` for
                                                         `LSZienkiewiczZhuEstimator`, and `Kee` for `KellyErrorEstimator`.
                                                            - When multiple substituted or overload types need to be encoded, append the
                                                            shorthands in argument-order as one PascalCase tail. Examples:
                                                            `CMFEM_ZienkiewiczZhuEstimator_NewDiGfFes`,
                                                            `CMFEM_ElasticityIntegrator_NewPwcPwc`,
                                                            `CMFEM_BilinearForm_FormLinearSystemSm`,
                                                            `CMFEM_GridFunction_ProjectBdrCoefficientFc`, and
                                                            `CMFEM_ThresholdRefiner_NewKee`.
                                                            - New wrapper types should expose the standard lifetime surface whenever the
                                                            underlying MFEM type is constructible:
                                                            `Construct`, `New`, `Copy`, `NewCopy`, `Destroy`, and `Delete`.
                                                            - `Construct` and `Copy` return the storage struct by value.
                                                            - `New` and `NewCopy` return pointers backed by actual C++ `new` allocations.
                                                            - `Destroy` manually invokes the C++ destructor on the storage object.
                                                            - `Delete` must be paired with `New`/`NewCopy` and call C++ `delete`.
                                                            - Because the public storage is only `char data[N]`, by-value construction in
                                                               the `.cpp` side must use an `alignas(mfem::<CppType>)` local before
                                                               placement-new into that storage.
                                                               - For the same alignment reason, C code that keeps one of these storage objects
                                                               by value on the stack should declare it with `_Alignas(max_align_t)`, as the
                                                               existing C examples do.
- Keep the wrapper surface minimal. Only wrap the MFEM APIs actually needed by
   the current C examples unless a broader wrapper expansion is explicitly
   requested.
- When adding a new wrapper, preserve the mirrored directory layout in
  `c_api/cmfem.h` include paths too. The aggregate header should include only
  the fundamental `c_api` wrapper surface, e.g. `fem/grid_function.h`.
- `c_api/cmfem.h` is a public package surface. It must never include, directly
  or transitively, headers from `examples/`.
- C examples that need example-specific bridge APIs must include the relevant
  `examples/c/adapters/exNN/...` headers directly instead of relying on
  `cmfem.h` to export them.
                                                               - Runtime path lookup for installed data is not hardcoded in public MFEM
                                                               headers. The canonical implementation lives in `config/runtime_paths.c` and
                                                               `config/runtime_paths.h`, with a C++ adapter in `config/runtime_paths.hpp`.
                                                               - C wrappers and C examples should call the runtime-path APIs instead of
                                                               embedding repo-relative paths.
                                                               - When wrapping MFEM's free iterative-solver helpers, normalize any odd
                                                               tolerance conventions to match the class-solver APIs that examples usually
                                                                  follow. In particular, the free `CG(...)`, `PCG(...)`, `GMRES(...)`, and
                                                                  `MINRES(...)` entry points use squared tolerance conventions relative to the
                                                                  corresponding solver-class APIs; the `c_api` wrappers should accept the
                                                                  ordinary relative/absolute tolerances that `CGSolver`, `GMRESSolver`, and
                                                                  `MINRESSolver` callers expect and square them internally.
                                                                     - Recent mixed/block wrappers also use these shorthands:
                                                                     `Bop` for `BlockOperator`, `Bdp` for
                                                                        `BlockDiagonalPreconditioner`, `Cgs` for `CGSolver`, `Rop` for
                                                                              `RAPOperator`, `RtTrace` for `RT_Trace_FECollection`, and `Dmg` for
                                                                                    the focused diffusion-multigrid helper surface.
                                                                                    - When an example only needs a narrow nonlinear assembly pattern, prefer a
                                                                                    focused helper free function in `c_api/` over exposing a broad generic
                                                                                    callback surface. The `ex40` support follows that pattern with
                                                                                    `CMFEM_AssembleIsomorphismLinearFormRtGf` and
                                                                                    `CMFEM_AssembleDIsomorphismMassMatrixRtGf`, instead of adding general
                                                                                    `ElementTransformation`-level coefficient callbacks to the public C ABI.
                                                                                    - The same rule applies to shared C++ example helper code that is not part of
                                                                                    MFEM's public API. If a C port only needs a narrow slice of helper logic,
                                                                                    wrap it through a focused `c_api/` bridge instead of reimplementing dense or
                                                                                    LAPACK-side algorithms in C. `ex33` follows that pattern with
                                                                                    `CMFEM_ComputePartialFractionApproximation`. If such a helper returns
                                                                                    heap-allocated buffers to C, expose an explicit matching free function such as
                                                                                    `CMFEM_FreeDoubles`.
                                                                                    - The same focused-bridge rule also applies to example-specific mesh builders
                                                                                    and postprocessing helpers when the current public C surface does not expose
                                                                                    the required low-level MFEM mutators cleanly. `ex27` follows that pattern
                                                                                    with `CMFEM_NewPeriodicTwoHoleMesh2d` and
                                                                                    `CMFEM_IntegrateBoundaryConditionGfAi`, while keeping the actual PDE assembly
                                                                                       and solve path in the C example.
                                                                                       - The same pattern is the right default for example-specific nonlinear
                                                                                          operators and Newton/preconditioner stacks. If a C port needs custom C++
                                                                                          classes like `HyperelasticOperator`, `ReducedSystemOperator`, or a block
                                                                                          `RubberOperator`, keep them behind a narrow helper context in `c_api/`
                                                                                          instead of wrapping broad generic `NonlinearForm`, `BlockNonlinearForm`,
                                                                                          `NewtonSolver`, or monitor APIs. `ex10` and `ex19` follow that pattern.
                                                                                          - When a focused helper owns an MFEM space hierarchy, match the ownership
                                                                                          expectations of the corresponding C++ example instead of adding extra
                                                                                          deletes on the C side. `ex26` follows the serial example's pattern:
                                                                                          the hierarchy helper owns the mesh/fespace chain, so the C example must
                                                                                          not separately delete the original mesh after deleting the multigrid
                                                                                          context.
                                                                                          - For saddle-point and other block-structured C example ports, prefer
                                                                                          explicit block offsets with flat `CMFEM_Vector` storage unless the example
                                                                                          truly needs `BlockVector` alias semantics. That keeps the public C surface
                                                                                          smaller and avoids wrapping MFEM container classes that are only being used
                                                                                          as temporary views.
                                                                                          - For MFEM ODE wrappers, keep the public C type as an adapter/handle around the
                                                                                          MFEM solver/operator object rather than moving example logic into C++.
                                                                                          The first-order and second-order operator wrappers in `c_api/` are callback
                                                                                          adapters; example-specific PDE logic should stay in the C example as a C
context struct plus callbacks.
- Two small generic wrappers added for the elasticity ports are easy to miss:
   `CMFEM_Mesh_GetNodesGf` is the C-side equivalent of `Mesh::GetNodes(
      GridFunction&)`
      and is the right way to recover reference coordinates on meshes that do not
         already own nodal grid functions; and
`CMFEM_FiniteElementSpace_NewMeshH1VDimByVdim` must be used when a port depends
on `Ordering::byVDIM` parity with the C++ example, such as `ex19`.

## Core Architecture Notes

- `core/` is the real MFEM public header root now. The old repo-root
   forwarding `mfem.hpp` was removed.
   - Code that still writes `#include "mfem.hpp"` relies on the build system
   adding `core/` to the include path so that `core/mfem.hpp` is found as the
   include-root `mfem.hpp`.
   - `config/` is a separate include root. Build systems must add both `core/`
   and `config/` include directories.
   - `core/mfem.hpp` includes `"config.hpp"`, not `"config/config.hpp"`. That only
   works when `config/` is an include root.
   - `config/config.hpp` supports out-of-tree configuration via
   `MFEM_CONFIG_FILE`. The current Zig build generates `config/_config.hpp` and
   points `MFEM_CONFIG_FILE` at it. Do not reintroduce the removed
   `core/config` symlink approach.
   - The umbrella include order in `core/mfem.hpp` is deliberate:
   `general -> linalg -> mesh -> fem`.
   - When possible, preserve that layering. Avoid creating new reverse
   dependencies from `general` into higher-level subsystems, or from `linalg`
   into `mesh`/`fem`.
   - The subsystem umbrellas are:
   `core/linalg/linalg.hpp`, `core/mesh/mesh_headers.hpp`, and `core/fem/fem.hpp`.
   If a new public header belongs to one of those surfaces, add it to the
   relevant umbrella instead of relying only on direct includes from examples.
   - `general/` is not just miscellaneous utilities. `device.hpp`,
   `forall.hpp`, `mem_manager.hpp`, `mem_alloc.hpp`, `globals.hpp`, and
   `error.hpp` form a cross-cutting substrate used throughout the library.
   - `mfem::Device` in `core/general/device.hpp` is effectively a singleton-style
   global configuration point. Many MFEM objects assume device configuration
   happens once early and is not later reconfigured.
   - OpenMP in `Device` is still a host backend. It is not a substitute for the
   MPI/distributed parts of MFEM.
   - Serial and parallel APIs are interleaved in the same source trees. The
   parallel surface is usually gated by `MFEM_USE_MPI` and often uses `p*`
   naming (`pmesh`, `pfespace`, `pgridfunc`, `plinearform`, `pbilinearform`,
   and related headers).
      - Optional TPL and device integrations also live in the same trees and are
      selected mostly by `MFEM_USE_*` macros rather than by separate top-level
      source roots.
      - `core/linalg/handle.hpp` defines `OperatorHandle`, and `OperatorPtr` is only
      a typedef of that class. Treat it like an owning/non-owning handle with
      `own_oper` semantics, not like a raw pointer wrapper.
   - That `OperatorPtr` detail matters in the C wrapper: copy operations are
      shallow unless you explicitly clone the underlying operator yourself.
   - `core/linalg/linalg.hpp` is where the serial/parallel split becomes obvious:
            base linear algebra is always present, while hypre/PETSc/SLEPc/MUMPS/etc. are
            added only under `MFEM_USE_MPI` or other feature macros.
         - `core/mesh/mesh_headers.hpp` and `core/fem/fem.hpp` follow the same pattern:
            serial headers are always present, while MPI-only or optional data-collection
               layers are conditionally included.
               - `core/fem/fem.hpp` unconditionally includes the serial finite-element stack
               and also `dfem/doperator.hpp`; the parallel finite-element classes are added
later under `MFEM_USE_MPI`.
- When porting build systems, it is usually better to turn feature macros on or
   off and adjust source lists accordingly than to edit public headers to hide
   unavailable features.
