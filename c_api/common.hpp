#ifndef CMFEM_COMMON_HPP
#define CMFEM_COMMON_HPP

#include "cmfem.h"

#include "config/runtime_paths.h"
#include "core/mfem.hpp"

#include <new>
#include <utility>

namespace cmfem
{

template <typename CPP, typename C>
inline CPP *As(C *value)
{
   return reinterpret_cast<CPP *>(value);
}

template <typename CPP, typename C>
inline const CPP *As(const C *value)
{
   return reinterpret_cast<const CPP *>(value);
}

#define CMFEM_ASSERT_TYPE(CType, CPPType) static_assert(sizeof(CType) == sizeof(CPPType))

inline mfem::Array<int> &ArrayIntRef(CMFEM_ArrayInt *array)
{
   return *As<mfem::Array<int> >(array);
}

inline const mfem::Array<int> &ArrayIntRef(const CMFEM_ArrayInt *array)
{
   return *As<const mfem::Array<int> >(array);
}

inline mfem::Vector &VectorRef(CMFEM_Vector *vector)
{
   return *As<mfem::Vector>(vector);
}

inline const mfem::Vector &VectorRef(const CMFEM_Vector *vector)
{
   return *As<const mfem::Vector>(vector);
}

inline mfem::OperatorPtr &OperatorPtrRef(CMFEM_OperatorPtr *op)
{
   return *As<mfem::OperatorPtr>(op);
}

inline const mfem::OperatorPtr &OperatorPtrRef(const CMFEM_OperatorPtr *op)
{
   return *As<const mfem::OperatorPtr>(op);
}

inline mfem::SparseMatrix &SparseMatrixFromOperator(CMFEM_OperatorPtr *op)
{
   return dynamic_cast<mfem::SparseMatrix &>(*OperatorPtrRef(op));
}

inline const mfem::SparseMatrix &SparseMatrixFromOperator(
   const CMFEM_OperatorPtr *op)
{
   return dynamic_cast<const mfem::SparseMatrix &>(*OperatorPtrRef(op));
}

} // namespace cmfem

#endif
