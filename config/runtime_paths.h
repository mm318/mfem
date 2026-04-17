#ifndef MFEM_RUNTIME_PATHS_C_H
#define MFEM_RUNTIME_PATHS_C_H

#include <stddef.h>

#if defined(__cplusplus)
extern "C" {
#endif

const char *MFEM_Runtime_ExamplesDataDir(void);
const char *MFEM_Runtime_ExamplesDataPath(const char *leaf);
size_t MFEM_Runtime_ExamplesDataPrefixLength(void);

const char *MFEM_Runtime_MiniappsDataDir(void);
const char *MFEM_Runtime_MiniappsDataPath(const char *leaf);
size_t MFEM_Runtime_MiniappsDataPrefixLength(void);

const char *MFEM_Runtime_TestsDataDir(void);
const char *MFEM_Runtime_TestsDataPath(const char *leaf);
size_t MFEM_Runtime_TestsDataPrefixLength(void);

const char *MFEM_Runtime_TestsMiniappsDir(void);
const char *MFEM_Runtime_TestsMiniappsPath(const char *leaf);
size_t MFEM_Runtime_TestsMiniappsPrefixLength(void);

#if defined(__cplusplus)
}
#endif

#endif
