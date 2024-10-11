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

    // LSH Stuff
    uint64_t T_;
    uint64_t K_;
    vector<vector<list<ListPoint *>>> hash_tables_; // hash_table [ <list_of_list_points>, ... ]; length of this should be t.
    vector<vector<pair<int, int>>> hyperplanes_;  // hyperplanes[t][k] = {0 ... N-1}^2, two random coordinates.
    void InitLSH();
    uint64_t num_vec_in_tables_ = 0;
    uint64_t max_num_vec_in_tables = 0;
    uint64_t num_hash_values_per_t_;
    uint64_t FindLSH(const ListPoint *lp, uint64_t t); // 0 ... 2^K-1, which bucket does *lp belong to in table t?
    void AddToLSHTable(ListPoint *lp);
    void RemoveFromLSHTable(const ListPoint *lp);
    ListPoint *min_vector_ = nullptr;

    int UpdateList(ListPoint* p);
    KleinSampler* sampler_;
    long n_; // Rank
    long m_; // Dim
    int64 best_sqr_norm_;
    int64 goal_sqr_norm_;
    queue<ListPoint*> queue_;
    // Statistics
    long iterations_;
    long collisions_;
    bool verbose_;
    long iteration_at_min_;
    // Adi:
    vector<tuple<long, long, long>> search_stats_1_; // iteration, table size, num_solns
    vector<tuple<long, long, long>> search_stats_2_; // iteration, table size, num_solns
    // END 
};

#endif
