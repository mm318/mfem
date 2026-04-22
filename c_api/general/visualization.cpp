#include "common.hpp"

extern "C" {

   void CMFEM_SendSolutionToGLVis(const CMFEM_Mesh *mesh,
                                  const CMFEM_GridFunction *grid_function,
                                  const char *host,
                                  int port)
   {
      mfem::socketstream sock(host, port);
      sock.precision(8);
      sock << "solution\n"
           << *cmfem::As<const mfem::Mesh>(mesh)
           << *cmfem::As<const mfem::GridFunction>(grid_function)
           << std::flush;
   }

   void CMFEM_SendMeshToGLVis(const CMFEM_Mesh *mesh,
                              const char *host,
                              int port)
   {
      mfem::socketstream sock(host, port);
      sock.precision(8);
      sock << "mesh\n"
           << *cmfem::As<const mfem::Mesh>(mesh)
           << std::flush;
   }

} // extern "C"
