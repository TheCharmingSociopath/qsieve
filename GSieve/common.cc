#include "common.h"

// NTL_vector_impl(int64, vec_int64)
// NTL_io_vector_impl(int64, vec_int64)
// NTL_vector_impl(vec_int64, vec_vec_int64)
// NTL_matrix_impl(int64, vec_int64, vec_vec_int64 ,mat_int64)
// NTL_vector_impl(vec_double, vec_vec_double)
// NTL_matrix_impl(double, vec_double, vec_vec_double ,mat_double)

ListPoint* NewListPoint(long dims) {
  ListPoint* p = new ListPoint;
  p->norm = 0;
  p->v.SetLength(dims);
  for(int i = 0; i < dims; ++i) {
    p->v[i] = 0;
  }
  return p;
}

void DeleteListPoint(ListPoint* p) {
  delete p;
}

void VecZZToListPoint(const NTL::vec_ZZ &v, ListPoint* p) {
  long dims = v.length();
  p->v.SetLength(dims);
  p->norm = 0;
  for(int i = 0; i < dims; ++i) {
    p->v[i] = to_long(v[i]);
    p->norm += p->v[i] * p->v[i];
  }
}

bool CheckReduceCount(const ListPoint* p1, const ListPoint* p2) {
  long dims = p1->v.length();
  int64 dot = 0;
  for (long i = 0; i < dims; ++i) {
    dot += p1->v[i] * p2->v[i];
  }
  if (abs(2 * dot) <= p2->norm)
    return false;
  return true;
}


bool Reduce(ListPoint* p1, const ListPoint* p2) {
  long dims = p1->v.length();
  int64 dot = 0;
  for (long i = 0; i < dims; ++i) {
    dot += p1->v[i] * p2->v[i];
  }
  if (abs(2 * dot) <= p2->norm)
    return false;
  long q = round((double) dot / p2->norm);
  for (long i = 0; i < dims; ++i) {
    p1->v[i] -= q * p2->v[i];
  }
  p1->norm = p1->norm + q * q * p2->norm -
      2 * q * dot;
  return true;
}

void MatInt64FromMatZZ(mat_int64& A, NTL::mat_ZZ B) {
  long rows = B.NumRows();
  long cols = B.NumCols();
  A.SetDims(rows, cols);
  for (long row = 0; row < rows; ++row) {
    for (long col = 0; col < cols; ++col) {
      A[row][col] = to_long(B[row][col]);
    }
  }
}

void MatDoubleFromMatRR(mat_double& A, NTL::mat_RR B) {
  long rows = B.NumRows();
  long cols = B.NumCols();
  A.SetDims(rows, cols);
  for (long row = 0; row < rows; ++row) {
    for (long col = 0; col < cols; ++col) {
      A[row][col] = to_double(B[row][col]);
    }
  }
}

