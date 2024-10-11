#ifndef __COMMON__
#define __COMMON__

#include <NTL/ZZ.h>
#include <NTL/vec_double.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>

typedef int64_t int64;

// Adi
// NTL_CLIENT
// NTL_vector_decl(int64, vec_int64)
// NTL_io_vector_decl(int64, vec_int64)
// NTL_vector_decl(vec_int64, vec_vec_int64)
// NTL_matrix_decl(int64, vec_int64, vec_vec_int64, mat_int64)
// NTL_vector_decl(vec_double, vec_vec_double)
// NTL_matrix_decl(double, vec_double, vec_vec_double, mat_double)

typedef NTL::Vec<int64> vec_int64;
typedef NTL::Vec<vec_int64> vec_vec_64;
typedef NTL::Vec<NTL::vec_double> vec_vec_double;
typedef NTL::Mat<int64> mat_int64;
typedef NTL::Mat<double> mat_double;

// Adi

// It is efficient to keep the square norm of each point
struct ListPoint {
  vec_int64 v;
  int64 norm;
};

ListPoint* NewListPoint(long dims);
void DeleteListPoint(ListPoint* p);
void VecZZToListPoint(const NTL::vec_ZZ &v, ListPoint* p);

void MatInt64FromMatZZ(mat_int64& A, NTL::mat_ZZ B);
void MatDoubleFromMatRR(mat_double& A, NTL::mat_RR B);

bool Reduce(ListPoint* p1, const ListPoint* p2);
bool CheckReduceCount(const ListPoint* p1, const ListPoint* p2);

#endif //__COMMON__
