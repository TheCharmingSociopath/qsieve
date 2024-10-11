// Sampling follows the work in:
// "Trapdoors for Hard Lattices and New Cryptographic Constructions"
// Craig Gentry, Chris Peikert, Vinod Vaikuntanathan
// See there for explaination of the variables s, s_prime, t etc...

#ifndef __SAMPLER__
#define __SAMPLER__

#include "common.h"
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>

NTL_CLIENT

class KleinSampler {
  public:
    void Init(const mat_ZZ &B);
    ListPoint* Sample();
  private:
    long n_;
    long m_;
    double t_;
    mat_int64 B_;
    mat_double mu_;
    vec_double coef_;
    vec_double s_prime_square_;
    long SampleZ(double c, double s_square);
};

#endif //__SAMPLER__
