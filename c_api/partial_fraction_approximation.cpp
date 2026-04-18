#include "common.hpp"

#include "examples/cpp/ex33.hpp"

#include <cstdlib>

extern "C" {

   void CMFEM_ComputePartialFractionApproximation(double *alpha,
                                                  int *size,
                                                  double **coeffs,
                                                  double **poles)
   {
      mfem::real_t alpha_value = static_cast<mfem::real_t>(*alpha);
      mfem::Array<mfem::real_t> coeffs_array;
      mfem::Array<mfem::real_t> poles_array;

      ComputePartialFractionApproximation(alpha_value, coeffs_array, poles_array);

      *alpha = static_cast<double>(alpha_value);
      *size = coeffs_array.Size();
      *coeffs = nullptr;
      *poles = nullptr;

      if (*size == 0)
      {
         return;
      }

      *coeffs = static_cast<double *>(std::malloc(sizeof(double) * (*size)));
      *poles = static_cast<double *>(std::malloc(sizeof(double) * (*size)));
      MFEM_VERIFY(*coeffs != nullptr && *poles != nullptr,
                  "Unable to allocate partial fraction buffers.");

      for (int i = 0; i < *size; i++)
      {
         (*coeffs)[i] = static_cast<double>(coeffs_array[i]);
         (*poles)[i] = static_cast<double>(poles_array[i]);
      }
   }

   void CMFEM_FreeDoubles(double *values)
   {
      std::free(values);
   }

} // extern "C"
