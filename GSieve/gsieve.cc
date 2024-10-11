#include "gsieve.h"
// Adi
#include <fstream>

// Cleans up used memory
void GSieve::CleanUp() {
  list<ListPoint*>::iterator lp_it;
  for (lp_it = list_.begin(); lp_it != list_.end(); ++lp_it)
    DeleteListPoint(*lp_it);
  list_.clear();
  while (!queue_.empty()) {
    DeleteListPoint(queue_.front());
    queue_.pop();
  }
}

// Initializes the reducer.
void GSieve::Init(const mat_ZZ &B, KleinSampler* sampler) {
  n_ = B.NumRows();
  m_ = B.NumCols();
  sampler_ = sampler;
  collisions_ = 0;
  iterations_ = 0;
  // Clean up list and queue
  CleanUp();
  // Here we should do LLL,BKZ reduction if not already done
  
  sampler_->Init(B);
  best_sqr_norm_ = to_long(B[0]*B[0]);
  ListPoint* p;
  long current_norm;
  for (int i = 0; i < n_; ++i) {
    p = NewListPoint(m_);
    VecZZToListPoint(B[i], p);
    current_norm = UpdateList(p);
    if (current_norm < best_sqr_norm_)
      best_sqr_norm_ = current_norm;
  }
  max_list_size_ = list_.size();
  goal_sqr_norm_ = 0;
  verbose_ = 0;
}

void GSieve::SetGoalSqrNorm(long norm) {
  goal_sqr_norm_ = norm;
}

void GSieve::SetVerbose(bool verbose) {
  verbose_ = verbose;
}

// Reduces recuirsively the point with all the points with smaller norm
// Adds the point to the list if we don't have colission
// and puts to the queue all the points with bigger norm
// that can be reduced with it.
int GSieve::UpdateList(ListPoint* p) {
  list<ListPoint*>::iterator lp_it, tmp_lp_it;
  ListPoint* lp;
  long lp_idx_counter = 0; // Adi
  bool needs_reduction = true;
  // Reduce the new lattice point
  // with the vectors with smaller norm
  while(needs_reduction) {
    needs_reduction = false;
    // Adi: Collect search_stats_1_
    //
    //
    long num_solns = 0;
    lp_idx_counter = 0;
    for (lp_it = list_.begin(); lp_it != list_.end(); ++lp_it) {
      lp = *lp_it;
      if (lp->norm > p->norm)
        break;
      ++lp_idx_counter;
      if (CheckReduceCount(p, lp)) {
              ++num_solns;
      }
    }
    search_stats_1_.push_back(make_tuple(iterations_, list_.size(), lp_idx_counter, num_solns));
    // Adi: END
    //
    for (lp_it = list_.begin(); lp_it != list_.end(); ++lp_it) {
      lp = *lp_it;
      if (lp->norm > p->norm)
        break;
      //If there is one reduction the vector should re-pass the list
      if (Reduce(p, lp)) {
        needs_reduction = true;
      }
    }
  }
  // We got a collision
  if (p->norm == 0) {
    DeleteListPoint(p);
    return 0;
  }
  // lp_it shows to the first point with bigger norm
  // this is where we will insert the new point
  list_.insert(lp_it, p);
  // Let's reduce now the vectors with bigger norm
  int num_solns2 = 0;

  tmp_lp_it = lp_it;
  while(tmp_lp_it != list_.end()) {
    if (CheckReduceCount(lp, p)) {
      ++num_solns2;
    }
    ++tmp_lp_it;
  }

  while (lp_it != list_.end()) {
    tmp_lp_it = lp_it;
    lp = *lp_it;
    ++lp_it;
    if (Reduce(lp, p)) {
      list_.erase(tmp_lp_it);
      queue_.push(lp);
    }
  }
  search_stats_2_.push_back(make_tuple(iterations_, list_.size(), lp_idx_counter, num_solns2));
  return p->norm;
}

bool GSieve::Start() {
  ListPoint* current_point;
  int64 current_norm;
  
  // Loop till you find a short enough vector,
  // or enough collisions.
  // The 0.1 * max_list_size_ + 200 is much more
  // than needed in general so we are almost sure
  // we have found the shortest vector.
   while ((best_sqr_norm_ > goal_sqr_norm_) && 
          (collisions_ < 0.2 * max_list_size_ + 200)) {
          // (collisions_ < 500)) { // Adi
    iterations_++;
    max_list_size_ = max(max_list_size_, long(list_.size()));
    // Get next point for reduction
    if (queue_.empty()) {
      current_point = sampler_->Sample();
    } else {
      current_point = queue_.front();
      queue_.pop();
    }
    // Reduce it and put the vectors of the list that
    // need reduction to the queue.
    current_norm = UpdateList(current_point);
    if (current_norm == 0) {
      collisions_++;
    }
    if (current_norm > 0 && current_norm < best_sqr_norm_) {
      best_sqr_norm_ = current_norm;
      iteration_at_min_ = iterations_;
    }
  }
  if (verbose_) {
    vec_int64 a;
    cout << std::endl;
    cout << "Dimension = " << n_ << "\n";
    cout << "Maximum list size = " << max_list_size_ << "\n";
    cout << "Number of iterations = " << iterations_ << "\n";
    cout << "Number of collisions = " << collisions_ << "\n";
    cout << "Iteration at which min found = " << iteration_at_min_ << std::endl;
  }
  cout << "Best vector = " << list_.front()->v << "\n";
  cout << "Square norm = " << best_sqr_norm_ << "\n";
  cout << std::endl;
  int tttt;
  cin >> tttt;
  std::ofstream outfile1, outfile2, outfile3;

  cout << "Search 1 Stats: iterations_, list size, lp_it idx, num_solns\n";
  for(int j = 0; j < search_stats_1_.size(); ++j) {
      cout << get<0>(search_stats_1_[j]) << " " << get<1>(search_stats_1_[j]) << " " << get<2>(search_stats_1_[j]) << " " << get<3>(search_stats_1_[j]) << endl;
  }
  
  cout << "Search 2 Stats: iterations_, list size, lp_it idx, num_solns\n";
  for(int j = 0; j < search_stats_2_.size(); ++j) {
      cout << get<0>(search_stats_2_[j]) << " " << get<1>(search_stats_2_[j]) << " " << get<2>(search_stats_2_[j]) << " " << get<3>(search_stats_2_[j])<< endl;
  }

  if (best_sqr_norm_ > goal_sqr_norm_)
    return false;
  return true;
}
