#if defined(__linux__)
#define _POSIX_C_SOURCE 200809L
#endif

#include "runtime_paths.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

typedef struct MFEM_RuntimePathCacheEntry
{
   char *key;
   char *value;
   struct MFEM_RuntimePathCacheEntry *next;
} MFEM_RuntimePathCacheEntry;

static char *mfem_runtime_examples_dir;
static char *mfem_runtime_miniapps_dir;
static char *mfem_runtime_tests_dir;
static char *mfem_runtime_tests_miniapps_dir;
static char *mfem_runtime_install_root_dir;
static MFEM_RuntimePathCacheEntry *mfem_runtime_examples_paths;
static MFEM_RuntimePathCacheEntry *mfem_runtime_miniapps_paths;
static MFEM_RuntimePathCacheEntry *mfem_runtime_tests_paths;
static MFEM_RuntimePathCacheEntry *mfem_runtime_tests_miniapps_paths;

static void MFEM_RuntimeFail(const char *message)
{
   fprintf(stderr, "MFEM runtime path error: %s\n", message);
   abort();
}

static char *MFEM_RuntimeStrdup(const char *text)
{
   const size_t len = strlen(text);
   char *copy = (char *)malloc(len + 1);
   if (!copy)
   {
      MFEM_RuntimeFail("out of memory");
   }
   memcpy(copy, text, len + 1);
   return copy;
}

static int MFEM_RuntimeIsSeparator(char ch)
{
   return ch == '/' || ch == '\\';
}

static char *MFEM_RuntimeJoin(const char *lhs, const char *rhs)
{
   const size_t lhs_len = strlen(lhs);
   const size_t rhs_len = strlen(rhs);
   const int need_sep =
      (lhs_len > 0 && rhs_len > 0 && !MFEM_RuntimeIsSeparator(lhs[lhs_len - 1]));
   char *joined = (char *)malloc(lhs_len + rhs_len + (size_t)need_sep + 1);
   if (!joined)
   {
      MFEM_RuntimeFail("out of memory");
   }

   memcpy(joined, lhs, lhs_len);
   if (need_sep)
   {
      joined[lhs_len] = '/';
   }
   memcpy(joined + lhs_len + (size_t)need_sep, rhs, rhs_len + 1);
   return joined;
}

static char *MFEM_RuntimeNormalize(const char *path)
{
   const size_t path_len = strlen(path);
   char *scratch = MFEM_RuntimeStrdup(path);
   char **components = (char **)malloc((path_len + 1) * sizeof(char *));
   size_t *component_lengths = (size_t *)malloc((path_len + 1) * sizeof(size_t));
   size_t component_count = 0;
   size_t prefix_len = 0;
   int absolute = 0;
   size_t i = 0;
   char *result;
   size_t result_len = 0;
   size_t write_at = 0;

   if (!components || !component_lengths)
   {
      MFEM_RuntimeFail("out of memory");
   }

   if (path_len >= 2 && scratch[1] == ':')
   {
      prefix_len = 2;
      i = 2;
      if (MFEM_RuntimeIsSeparator(scratch[i]))
      {
         absolute = 1;
         while (MFEM_RuntimeIsSeparator(scratch[i]))
         {
            i++;
         }
      }
   }
   else if (MFEM_RuntimeIsSeparator(scratch[0]))
   {
      absolute = 1;
      while (MFEM_RuntimeIsSeparator(scratch[i]))
      {
         i++;
      }
   }

   while (scratch[i] != '\0')
   {
      char *component = scratch + i;
      size_t component_len = 0;
      while (scratch[i] != '\0' && !MFEM_RuntimeIsSeparator(scratch[i]))
      {
         i++;
         component_len++;
      }
      while (MFEM_RuntimeIsSeparator(scratch[i]))
      {
         scratch[i] = '\0';
         i++;
      }

      if (component_len == 0 ||
          (component_len == 1 && component[0] == '.'))
      {
         continue;
      }

      if (component_len == 2 && component[0] == '.' && component[1] == '.')
      {
         if (component_count > 0 &&
             !(component_lengths[component_count - 1] == 2 &&
               components[component_count - 1][0] == '.' &&
               components[component_count - 1][1] == '.'))
         {
            component_count--;
            continue;
         }

         if (!absolute)
         {
            components[component_count] = component;
            component_lengths[component_count] = component_len;
            component_count++;
         }
         continue;
      }

      components[component_count] = component;
      component_lengths[component_count] = component_len;
      component_count++;
   }

   result_len += prefix_len;
   if (absolute)
   {
      result_len++;
   }
   if (component_count == 0 && !absolute && prefix_len == 0)
   {
      result_len = 1;
   }
   else
   {
      for (i = 0; i < component_count; i++)
      {
         if (i > 0)
         {
            result_len++;
         }
         result_len += component_lengths[i];
      }
   }

   result = (char *)malloc(result_len + 1);
   if (!result)
   {
      MFEM_RuntimeFail("out of memory");
   }

   if (prefix_len > 0)
   {
      memcpy(result, scratch, prefix_len);
      write_at += prefix_len;
   }
   if (absolute)
   {
      result[write_at++] = '/';
   }
   if (component_count == 0 && !absolute && prefix_len == 0)
   {
      result[write_at++] = '.';
   }
   else
   {
      for (i = 0; i < component_count; i++)
      {
         if (i > 0)
         {
            result[write_at++] = '/';
         }
         memcpy(result + write_at, components[i], component_lengths[i]);
         write_at += component_lengths[i];
      }
   }
   result[write_at] = '\0';

   free(component_lengths);
   free(components);
   free(scratch);
   return result;
}

static char *MFEM_RuntimeExecutablePath(void)
{
#if defined(_WIN32)
   DWORD size = 32768;
   char *buffer = (char *)malloc((size_t)size + 1);
   DWORD actual;
   if (!buffer)
   {
      MFEM_RuntimeFail("out of memory");
   }
   actual = GetModuleFileNameA(NULL, buffer, size);
   if (actual == 0 || actual == size)
   {
      free(buffer);
      MFEM_RuntimeFail("unable to detect executable path on Windows");
   }
   buffer[actual] = '\0';
   {
      char *normalized = MFEM_RuntimeNormalize(buffer);
      free(buffer);
      return normalized;
   }
#elif defined(__APPLE__)
   uint32_t size = 0;
   char *buffer;
   if (_NSGetExecutablePath(NULL, &size) != -1)
   {
      MFEM_RuntimeFail("unexpected macOS executable path size result");
   }
   buffer = (char *)malloc((size_t)size + 1);
   if (!buffer)
   {
      MFEM_RuntimeFail("out of memory");
   }
   if (_NSGetExecutablePath(buffer, &size) != 0)
   {
      free(buffer);
      MFEM_RuntimeFail("unable to detect executable path on macOS");
   }
   {
      char *normalized = MFEM_RuntimeNormalize(buffer);
      free(buffer);
      return normalized;
   }
#elif defined(__linux__)
   size_t size = 4096;
   char *buffer = (char *)malloc(size);
   while (1)
   {
      ssize_t actual;
      if (!buffer)
      {
         MFEM_RuntimeFail("out of memory");
      }
      actual = readlink("/proc/self/exe", buffer, size - 1);
      if (actual < 0)
      {
         free(buffer);
         MFEM_RuntimeFail("unable to read /proc/self/exe");
      }
      if ((size_t)actual < size - 1)
      {
         buffer[actual] = '\0';
         {
            char *normalized = MFEM_RuntimeNormalize(buffer);
            free(buffer);
            return normalized;
         }
      }
      size *= 2;
      buffer = (char *)realloc(buffer, size);
   }
#else
   MFEM_RuntimeFail("unsupported platform for executable path lookup");
   return NULL;
#endif
}

static const char *MFEM_RuntimeExecutableDir(void)
{
   static char *dir;
   if (!dir)
   {
      char *path = MFEM_RuntimeExecutablePath();
      char *last_sep = path + strlen(path);
      while (last_sep > path && !MFEM_RuntimeIsSeparator(last_sep[-1]))
      {
         last_sep--;
      }
      if (last_sep == path)
      {
         dir = MFEM_RuntimeStrdup(".");
      }
      else
      {
         while (last_sep > path && MFEM_RuntimeIsSeparator(last_sep[-1]))
         {
            last_sep--;
         }
         if (last_sep == path)
         {
            dir = MFEM_RuntimeStrdup("/");
         }
         else
         {
            const size_t len = (size_t)(last_sep - path);
            dir = (char *)malloc(len + 1);
            if (!dir)
            {
               MFEM_RuntimeFail("out of memory");
            }
            memcpy(dir, path, len);
            dir[len] = '\0';
         }
      }
      free(path);
   }
   return dir;
}

static const char *MFEM_RuntimeInstallRootDir(void)
{
   if (!mfem_runtime_install_root_dir)
   {
      const char *dir = MFEM_RuntimeExecutableDir();
      const char *last_match = NULL;
      const char *cursor = dir;

      while ((cursor = strchr(cursor, '/')) != NULL)
      {
         if (strncmp(cursor, "/bin", 4) == 0 &&
             (cursor[4] == '\0' || cursor[4] == '/'))
         {
            last_match = cursor;
         }
         cursor++;
      }

      if (!last_match)
      {
         MFEM_RuntimeFail("executable is not installed under zig-out/bin");
      }

      if (last_match == dir)
      {
         mfem_runtime_install_root_dir = MFEM_RuntimeStrdup("/");
      }
      else
      {
         const size_t len = (size_t)(last_match - dir);
         mfem_runtime_install_root_dir = (char *)malloc(len + 1);
         if (!mfem_runtime_install_root_dir)
         {
            MFEM_RuntimeFail("out of memory");
         }
         memcpy(mfem_runtime_install_root_dir, dir, len);
         mfem_runtime_install_root_dir[len] = '\0';
      }
   }
   return mfem_runtime_install_root_dir;
}

static const char *MFEM_RuntimeCachedDir(const char *install_relative_path,
                                         char **slot)
{
   if (!*slot)
   {
      char *joined =
         MFEM_RuntimeJoin(MFEM_RuntimeInstallRootDir(), install_relative_path);
      *slot = MFEM_RuntimeNormalize(joined);
      free(joined);
   }
   return *slot;
}

static const char *MFEM_RuntimeCachedPath(const char *relative_dir,
                                          const char *leaf,
                                          char **dir_slot,
                                          MFEM_RuntimePathCacheEntry **head)
{
   MFEM_RuntimePathCacheEntry *entry = *head;
   const char *dir = MFEM_RuntimeCachedDir(relative_dir, dir_slot);
   size_t key_len = strlen(leaf);
   char *key;

   while (entry)
   {
      if (strcmp(entry->key, leaf) == 0)
      {
         return entry->value;
      }
      entry = entry->next;
   }

   entry = (MFEM_RuntimePathCacheEntry *)malloc(sizeof(*entry));
   key = (char *)malloc(key_len + 1);
   if (!entry || !key)
   {
      MFEM_RuntimeFail("out of memory");
   }
   memcpy(key, leaf, key_len + 1);
   entry->key = key;
   {
      char *joined = MFEM_RuntimeJoin(dir, leaf);
      entry->value = MFEM_RuntimeNormalize(joined);
      free(joined);
   }
   entry->next = *head;
   *head = entry;
   return entry->value;
}

const char *MFEM_Runtime_ExamplesDataDir(void)
{
   return MFEM_RuntimeCachedDir("share/data", &mfem_runtime_examples_dir);
}

const char *MFEM_Runtime_ExamplesDataPath(const char *leaf)
{
   return MFEM_RuntimeCachedPath("share/data",
                                 leaf,
                                 &mfem_runtime_examples_dir,
                                 &mfem_runtime_examples_paths);
}

size_t MFEM_Runtime_ExamplesDataPrefixLength(void)
{
   return strlen(MFEM_Runtime_ExamplesDataDir()) + 1;
}

const char *MFEM_Runtime_MiniappsDataDir(void)
{
   return MFEM_RuntimeCachedDir("share/data", &mfem_runtime_miniapps_dir);
}

const char *MFEM_Runtime_MiniappsDataPath(const char *leaf)
{
   return MFEM_RuntimeCachedPath("share/data",
                                 leaf,
                                 &mfem_runtime_miniapps_dir,
                                 &mfem_runtime_miniapps_paths);
}

size_t MFEM_Runtime_MiniappsDataPrefixLength(void)
{
   return strlen(MFEM_Runtime_MiniappsDataDir()) + 1;
}

const char *MFEM_Runtime_TestsDataDir(void)
{
   return MFEM_RuntimeCachedDir("share/data", &mfem_runtime_tests_dir);
}

const char *MFEM_Runtime_TestsDataPath(const char *leaf)
{
   return MFEM_RuntimeCachedPath("share/data",
                                 leaf,
                                 &mfem_runtime_tests_dir,
                                 &mfem_runtime_tests_paths);
}

size_t MFEM_Runtime_TestsDataPrefixLength(void)
{
   return strlen(MFEM_Runtime_TestsDataDir()) + 1;
}

const char *MFEM_Runtime_TestsMiniappsDir(void)
{
   return MFEM_RuntimeCachedDir("share/miniapps",
                                &mfem_runtime_tests_miniapps_dir);
}

const char *MFEM_Runtime_TestsMiniappsPath(const char *leaf)
{
   return MFEM_RuntimeCachedPath("share/miniapps",
                                 leaf,
                                 &mfem_runtime_tests_miniapps_dir,
                                 &mfem_runtime_tests_miniapps_paths);
}

size_t MFEM_Runtime_TestsMiniappsPrefixLength(void)
{
   return strlen(MFEM_Runtime_TestsMiniappsDir()) + 1;
}
