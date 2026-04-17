#ifndef MFEM_RUNTIME_PATHS_HPP
#define MFEM_RUNTIME_PATHS_HPP

#include "runtime_paths.h"

#include <string>

namespace mfem
{
namespace test
{

inline const std::string &ExamplesDataDir()
{
   static const std::string dir = MFEM_Runtime_ExamplesDataDir();
   return dir;
}

inline const char *ExamplesDataPath(const char *leaf)
{
   return MFEM_Runtime_ExamplesDataPath(leaf);
}

inline std::string ExamplesDataPath(const std::string &leaf)
{
   return MFEM_Runtime_ExamplesDataPath(leaf.c_str());
}

inline std::size_t ExamplesDataPrefixLength()
{
   return MFEM_Runtime_ExamplesDataPrefixLength();
}

inline const std::string &MiniappsDataDir()
{
   static const std::string dir = MFEM_Runtime_MiniappsDataDir();
   return dir;
}

inline const char *MiniappsDataPath(const char *leaf)
{
   return MFEM_Runtime_MiniappsDataPath(leaf);
}

inline std::string MiniappsDataPath(const std::string &leaf)
{
   return MFEM_Runtime_MiniappsDataPath(leaf.c_str());
}

inline std::size_t MiniappsDataPrefixLength()
{
   return MFEM_Runtime_MiniappsDataPrefixLength();
}

inline const std::string &TestsDataDir()
{
   static const std::string dir = MFEM_Runtime_TestsDataDir();
   return dir;
}

inline const char *TestsDataPath(const char *leaf)
{
   return MFEM_Runtime_TestsDataPath(leaf);
}

inline std::string TestsDataPath(const std::string &leaf)
{
   return MFEM_Runtime_TestsDataPath(leaf.c_str());
}

inline std::size_t TestsDataPrefixLength()
{
   return MFEM_Runtime_TestsDataPrefixLength();
}

inline const std::string &TestsMiniappsDir()
{
   static const std::string dir = MFEM_Runtime_TestsMiniappsDir();
   return dir;
}

inline const char *TestsMiniappsPath(const char *leaf)
{
   return MFEM_Runtime_TestsMiniappsPath(leaf);
}

inline std::string TestsMiniappsPath(const std::string &leaf)
{
   return MFEM_Runtime_TestsMiniappsPath(leaf.c_str());
}

inline std::size_t TestsMiniappsPrefixLength()
{
   return MFEM_Runtime_TestsMiniappsPrefixLength();
}

} // namespace test
} // namespace mfem

#endif
