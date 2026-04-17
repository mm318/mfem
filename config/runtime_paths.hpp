// Runtime helpers for locating installed MFEM data relative to executables in
// zig-out/bin.
#ifndef MFEM_RUNTIME_PATHS_HPP
#define MFEM_RUNTIME_PATHS_HPP

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

namespace mfem
{
namespace test
{

namespace detail
{

[[noreturn]] inline void FailExecutablePathLookup(const char *message)
{
   std::fprintf(stderr, "MFEM runtime path error: %s\n", message);
   std::abort();
}

inline std::filesystem::path DetectExecutablePath()
{
#if defined(_WIN32)
   std::vector<char> buffer(32768);
   const DWORD size = GetModuleFileNameA(nullptr, buffer.data(),
                                         static_cast<DWORD>(buffer.size()));
   if (size == 0 || size == buffer.size())
   {
      FailExecutablePathLookup("unable to detect executable path on Windows");
   }
   return std::filesystem::path(std::string(buffer.data(), size));
#elif defined(__APPLE__)
   std::uint32_t size = 0;
   _NSGetExecutablePath(nullptr, &size);
   std::vector<char> buffer(size + 1);
   if (_NSGetExecutablePath(buffer.data(), &size) != 0)
   {
      FailExecutablePathLookup("unable to detect executable path on macOS");
   }
   return std::filesystem::path(buffer.data()).lexically_normal();
#elif defined(__linux__)
   std::vector<char> buffer(4096);
   while (true)
   {
      const ssize_t size = ::readlink("/proc/self/exe",
                                      buffer.data(),
                                      buffer.size() - 1);
      if (size < 0)
      {
         FailExecutablePathLookup("unable to read /proc/self/exe");
      }
      if (size < static_cast<ssize_t>(buffer.size() - 1))
      {
         buffer[size] = '\0';
         return std::filesystem::path(buffer.data()).lexically_normal();
      }
      buffer.resize(buffer.size() * 2);
   }
#else
   FailExecutablePathLookup("unsupported platform for executable path lookup");
#endif
}

inline const std::string &ExecutableDir()
{
   static const std::string dir =
      DetectExecutablePath().parent_path().lexically_normal().string();
   return dir;
}

inline const std::string &CachedInstalledDir(const char *relative_path)
{
   static std::mutex mutex;
   static std::unordered_map<std::string, std::string> cache;

   std::lock_guard<std::mutex> lock(mutex);
   auto [it, inserted] = cache.try_emplace(relative_path);
   if (inserted)
   {
      it->second =
         (std::filesystem::path(ExecutableDir()) / relative_path)
         .lexically_normal().string();
   }
   return it->second;
}

inline std::string InstalledPath(const char *relative_dir,
                                 const std::string &leaf)
{
   return (std::filesystem::path(CachedInstalledDir(relative_dir)) / leaf)
          .lexically_normal().string();
}

inline const char *CachedInstalledPath(const char *relative_dir,
                                       const char *leaf)
{
   static std::mutex mutex;
   static std::unordered_map<std::string, std::string> cache;

   const std::string key = std::string(relative_dir) + '\n' + leaf;

   std::lock_guard<std::mutex> lock(mutex);
   auto [it, inserted] = cache.try_emplace(key);
   if (inserted)
   {
      it->second = InstalledPath(relative_dir, leaf);
   }
   return it->second.c_str();
}

} // namespace detail

inline const std::string &ExamplesDataDir()
{
   return detail::CachedInstalledDir("../../../share/data");
}

inline const char *ExamplesDataPath(const char *leaf)
{
   return detail::CachedInstalledPath("../../../share/data", leaf);
}

inline std::string ExamplesDataPath(const std::string &leaf)
{
   return detail::InstalledPath("../../../share/data", leaf);
}

inline std::size_t ExamplesDataPrefixLength()
{
   return ExamplesDataDir().size() + 1;
}

inline const std::string &MiniappsDataDir()
{
   return detail::CachedInstalledDir("../../share/data");
}

inline const char *MiniappsDataPath(const char *leaf)
{
   return detail::CachedInstalledPath("../../share/data", leaf);
}

inline std::string MiniappsDataPath(const std::string &leaf)
{
   return detail::InstalledPath("../../share/data", leaf);
}

inline std::size_t MiniappsDataPrefixLength()
{
   return MiniappsDataDir().size() + 1;
}

inline const std::string &TestsDataDir()
{
   return detail::CachedInstalledDir("../../share/data");
}

inline const char *TestsDataPath(const char *leaf)
{
   return detail::CachedInstalledPath("../../share/data", leaf);
}

inline std::string TestsDataPath(const std::string &leaf)
{
   return detail::InstalledPath("../../share/data", leaf);
}

inline std::size_t TestsDataPrefixLength()
{
   return TestsDataDir().size() + 1;
}

inline const std::string &TestsMiniappsDir()
{
   return detail::CachedInstalledDir("../../share/miniapps");
}

inline const char *TestsMiniappsPath(const char *leaf)
{
   return detail::CachedInstalledPath("../../share/miniapps", leaf);
}

inline std::string TestsMiniappsPath(const std::string &leaf)
{
   return detail::InstalledPath("../../share/miniapps", leaf);
}

inline std::size_t TestsMiniappsPrefixLength()
{
   return TestsMiniappsDir().size() + 1;
}

} // namespace test
} // namespace mfem

#endif
