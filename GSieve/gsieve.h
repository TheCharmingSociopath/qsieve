#ifndef __GSIEVE__
#define __GSIEVE__

#include <iostream>
#include <list>
#include <queue>
#include <vector>
#include <tuple>

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include "common.h"
#include "sampler.h"

NTL_CLIENT

class GSieve {
  public:
    ~GSieve() {
      CleanUp();
    }
    void Init(const mat_ZZ& B, KleinSampler* sampler);
    void SetGoalSqrNorm(long norm);
    void SetVerbose(bool verbose);
    void CleanUp();
    bool Start();
  private:
    int UpdateList(ListPoint* p);
    KleinSampler* sampler_;
    long n_;
    long m_;
    int64 best_sqr_norm_;
    int64 goal_sqr_norm_;
    list<ListPoint*> list_;
    queue<ListPoint*> queue_;
    // Statistics
    long max_list_size_;
    long iterations_;
    long collisions_;
    bool verbose_;
    long iteration_at_min_;
    // Adi:
    vector<tuple<long, long, long, long>> search_stats_1_; // iteration, list size, lp_it idx, num_solns.  lp_it idx: list is sorted with norm...this is the point till which norm(v) <= current_norm. 
    vector<tuple<long, long, long, long>> search_stats_2_; // iteration, list size, lp_it idx, num_solns.  lp_it idx for the second search is the opposite, idx of first vector with norm bigger. Second search runs this point onwards. 
    // END 
};

#endif
