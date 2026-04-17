#include "common.hpp"

extern "C" {

   const char *CMFEM_ExamplesDataPath(const char *leaf)
   {
      return MFEM_Runtime_ExamplesDataPath(leaf);
   }

} // extern "C"
