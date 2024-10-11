#include "gsieve.h"
// Adi
#include <fstream>

uint64_t GSieve::FindLSH(const ListPoint *lp, uint64_t t)
{
  uint64_t ret = 0;

  for (int j = 0; j < K_; ++j)
  {
    ret = ret << 1;
    uint64_t ip = lp->v[hyperplanes_[t][j].first] - lp->v[hyperplanes_[t][j].second];
    if (ip > 0)
      ++ret;
  }

  // As h(-v) = -h(v), merge bucket u with 2^K_ - u - 1
  if (ret >= (1 << (K_ - 1)))
    ret = 2 * (1 << (K_ - 1)) - ret - 1;
  return ret;
}

void GSieve::InitLSH()
{
  hyperplanes_.resize(T_);
  hash_tables_.resize(T_);
  for (int i = 0; i < T_; ++i)
  {
    hyperplanes_[i].resize(K_);
    hash_tables_[i].resize(num_hash_values_per_t_);
    for (int j = 0; j < K_; ++j)
    {
      hyperplanes_[i][j].first = rand() % m_;
      hyperplanes_[i][j].second = rand() % m_;
      while (hyperplanes_[i][j].second == hyperplanes_[i][j].first)
        hyperplanes_[i][j].second = rand() % m_;
    }
  }
}

void GSieve::AddToLSHTable(ListPoint *lp)
{
  ++num_vec_in_tables_;
  max_num_vec_in_tables = max(max_num_vec_in_tables, num_vec_in_tables_);
  for (int i = 0; i < T_; ++i)
  {
    uint64_t hash = FindLSH(lp, i);
    // maintain sorted norm.
    auto lp_it = hash_tables_[i][hash].begin();
    while (lp_it != hash_tables_[i][hash].end() and (*lp_it)->norm < lp->norm)
    {
      ++lp_it;
    }
    hash_tables_[i][hash].insert(lp_it, lp);
  }
}

void GSieve::RemoveFromLSHTable(const ListPoint *lp)
{
  --num_vec_in_tables_;
  for (int i = 0; i < T_; ++i)
  {
    uint64_t hash = FindLSH(lp, i);
    for (auto lp_it = hash_tables_[i][hash].begin(); lp_it != hash_tables_[i][hash].end(); ++lp_it)
    {
      if ((*lp_it)->norm == lp->norm and (*lp_it)->v == lp->v)
      {
        hash_tables_[i][hash].erase(lp_it);
        break;
      }
    }
  }
}

// Cleans up used memory
void GSieve::CleanUp()
{
  for (int j = 0; j < hash_tables_.size(); ++j)
  {
    for (int i = 0; i < hash_tables_[j].size(); ++i)
    {
      auto lp_it = hash_tables_[j][i].begin();
      while (lp_it != hash_tables_[j][i].end())
      {
        auto temp_lp_ptr = *lp_it;
        ++lp_it;
        RemoveFromLSHTable(temp_lp_ptr);
        DeleteListPoint(temp_lp_ptr);
      }
      hash_tables_[j][i].clear();
    }
  }

  while (!queue_.empty())
  {
    DeleteListPoint(queue_.front());
    queue_.pop();
  }
}

// Initializes the reducer.
void GSieve::Init(const mat_ZZ &B, KleinSampler *sampler)
{
  n_ = B.NumRows(); // Rank
  m_ = B.NumCols(); // Dim
  sampler_ = sampler;
  collisions_ = 0;
  iterations_ = 0;
  // Clean up list and queue
  CleanUp();
  // Here we should do LLL,BKZ reduction if not already done

  sampler_->Init(B);
  T_ = ceil(powf(2, 0.1290F * m_));
  K_ = ceil(0.2209F * m_);
  num_hash_values_per_t_ = powf(2, K_ - 1);
  InitLSH();
  best_sqr_norm_ = to_long(B[0] * B[0]);
  ListPoint *p;
  long current_norm;
  for (int i = 0; i < n_; ++i)
  {
    p = NewListPoint(m_);
    VecZZToListPoint(B[i], p);
    current_norm = UpdateList(p);
    if (current_norm < best_sqr_norm_ or i == 0) {
      best_sqr_norm_ = current_norm;
      min_vector_ = p;
    }
  }
  goal_sqr_norm_ = 0;
  verbose_ = 0;
}

void GSieve::SetGoalSqrNorm(long norm)
{
  goal_sqr_norm_ = norm;
}

void GSieve::SetVerbose(bool verbose)
{
  verbose_ = verbose;
}

// Reduces recuirsively the point with all the points with smaller norm
// Adds the point to the list if we don't have colission
// and puts to the queue all the points with bigger norm
// that can be reduced with it.
int GSieve::UpdateList(ListPoint *p) // Now Update LSH table.
{
  bool needs_reduction = true;
  // Reduce the new lattice point
  // with the vectors with smaller norm
  while (needs_reduction)
  {
    needs_reduction = false;
    for (int t = 0; t < T_; ++t)
    {
      ListPoint *temp_lp_ptr;
      uint64_t hash = FindLSH(p, t);
      long num_solns = 0;
      bool soln_not_found = true;
      for (auto lp : hash_tables_[t][hash])
      {
        if (CheckReduceCount(p, lp))
        {
          ++num_solns;
          if (soln_not_found)
          {
            temp_lp_ptr = lp;
            soln_not_found = false;
          }
        }
      }
      search_stats_1_.push_back(make_tuple(iterations_, hash_tables_[t][hash].size(), num_solns));

      if (num_solns > 0)
      {
        Reduce(p, temp_lp_ptr);
        needs_reduction = true;
        break;
      }
    }
  }
  // We got a collision
  if (p->norm == 0)
  {
    DeleteListPoint(p);
    return 0;
  }

  // Let's reduce now the vectors with bigger norm
  for (int t = 0; t < T_; ++t)
  {
    uint64_t hash = FindLSH(p, t);
    queue<ListPoint *> temp_q;
    int num_solns2 = 0;
    auto lp_it = hash_tables_[t][hash].begin();

    while (lp_it != hash_tables_[t][hash].end())
    {
      auto prev_lp_ptr = *lp_it;
      ++lp_it;
      if (CheckReduceCount(prev_lp_ptr, p))
      {
        ++num_solns2;
        RemoveFromLSHTable(prev_lp_ptr);
        Reduce(prev_lp_ptr, p);
        queue_.push(prev_lp_ptr);
      }
    }
    search_stats_2_.push_back(make_tuple(iterations_, hash_tables_[t][hash].size(), num_solns2));
  }
  AddToLSHTable(p);
  return p->norm;
}

bool GSieve::Start()
{
  ListPoint *current_point;
  int64 current_norm;

  // Loop till you find a short enough vector,
  // or enough collisions.
  // The 0.1 * max_list_size_ + 200 is much more
  // than needed in general so we are almost sure
  // we have found the shortest vector.
  while ((best_sqr_norm_ > goal_sqr_norm_) &&
         (collisions_ < 0.1 * max_num_vec_in_tables + 200))
  //  (collisions_ < 500))
  { // Adi
    iterations_++;
    // Get next point for reduction
    if (queue_.empty())
    {
      current_point = sampler_->Sample();
    }
    else
    {
      current_point = queue_.front();
      queue_.pop();
    }
    // Reduce it and put the vectors of the list that
    // need reduction to the queue.
    current_norm = UpdateList(current_point);
    if (current_norm == 0)
    {
      collisions_++;
    }
    if (current_norm > 0 && current_norm < best_sqr_norm_)
    {
      min_vector_ = current_point;
      best_sqr_norm_ = current_norm;
      iteration_at_min_ = iterations_;
    }
  }
  if (verbose_)
  {
    vec_int64 a;
    cout << std::endl;
    cout << "Dimension = " << n_ << "\n";
    cout << "Maximum num vectors in tables = " << max_num_vec_in_tables << "\n";
    cout << "Number of iterations = " << iterations_ << "\n";
    cout << "Number of collisions = " << collisions_ << "\n";
    cout << "Iteration at which min found = " << iteration_at_min_ << std::endl;
  }

  cout << "Best vector = " << min_vector_->v << "\n";
  cout << "Square norm = " << min_vector_->norm << "\n";
  cout << std::endl;
  
  std::ofstream outfile1, outfile2, outfile3;

  cout << "Search 1 Stats: iterations_, table size, num_solns\n";
  for (int j = 0; j < search_stats_1_.size(); ++j)
  {
    cout << get<0>(search_stats_1_[j]) << " " << get<1>(search_stats_1_[j]) << " " << get<2>(search_stats_1_[j]) << endl;
  }

  cout << "Search 2 Stats: iterations_, table size, num_solns\n";
  for (int j = 0; j < search_stats_2_.size(); ++j)
  {
    cout << get<0>(search_stats_2_[j]) << " " << get<1>(search_stats_2_[j]) << " " << get<2>(search_stats_2_[j]) << endl;
  }

  if (best_sqr_norm_ > goal_sqr_norm_)
    return false;
  return true;
}
