#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_Vector, mfem::Vector);

} // namespace

extern "C" {

CMFEM_Vector CMFEM_Vector_Construct(void)
{
   alignas(mfem::Vector) CMFEM_Vector vector;
   new (cmfem::As<mfem::Vector>(&vector)) mfem::Vector();
   return vector;
}

CMFEM_Vector CMFEM_Vector_ConstructSize(int size)
{
   alignas(mfem::Vector) CMFEM_Vector vector;
   new (cmfem::As<mfem::Vector>(&vector)) mfem::Vector(size);
   return vector;
}

CMFEM_Vector *CMFEM_Vector_New(void)
{
   return reinterpret_cast<CMFEM_Vector *>(new mfem::Vector());
}

CMFEM_Vector *CMFEM_Vector_NewSize(int size)
{
   return reinterpret_cast<CMFEM_Vector *>(new mfem::Vector(size));
}

CMFEM_Vector CMFEM_Vector_Copy(const CMFEM_Vector *vector)
{
   alignas(mfem::Vector) CMFEM_Vector copy;
   new (cmfem::As<mfem::Vector>(&copy)) mfem::Vector(cmfem::VectorRef(vector));
   return copy;
}

CMFEM_Vector *CMFEM_Vector_NewCopy(const CMFEM_Vector *vector)
{
   return reinterpret_cast<CMFEM_Vector *>(new mfem::Vector(cmfem::VectorRef(vector)));
}

void CMFEM_Vector_Destroy(CMFEM_Vector *vector)
{
   cmfem::As<mfem::Vector>(vector)->~Vector();
}

void CMFEM_Vector_Delete(CMFEM_Vector *vector)
{
   delete cmfem::As<mfem::Vector>(vector);
}

void CMFEM_Vector_Assign(CMFEM_Vector *vector, double value)
{
   cmfem::VectorRef(vector) = static_cast<mfem::real_t>(value);
}

void CMFEM_Vector_Set(CMFEM_Vector *vector, int index, double value)
{
   cmfem::VectorRef(vector)(index) = static_cast<mfem::real_t>(value);
}

double CMFEM_Vector_Get(const CMFEM_Vector *vector, int index)
{
   return cmfem::VectorRef(vector)(index);
}

int CMFEM_Vector_Size(const CMFEM_Vector *vector)
{
   return cmfem::VectorRef(vector).Size();
}

} // extern "C"
