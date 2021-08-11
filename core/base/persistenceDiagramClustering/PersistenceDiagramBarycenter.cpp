#include <PersistenceDiagramBarycenter.h>

#include <numeric>

void ttk::PersistenceDiagramBarycenter::execute(
  std::vector<DiagramType> &intermediateDiagrams,
  DiagramType &barycenter,
  std::vector<std::vector<MatchingType>> &all_matchings) {

  Timer tm;

  printMsg("Computing Barycenter of " + std::to_string(numberOfInputs_)
           + " diagrams.");

  std::vector<DiagramType> data_min(numberOfInputs_);
  std::vector<DiagramType> data_sad(numberOfInputs_);
  std::vector<DiagramType> data_max(numberOfInputs_);

  std::vector<std::vector<int>> data_min_idx(numberOfInputs_);
  std::vector<std::vector<int>> data_sad_idx(numberOfInputs_);
  std::vector<std::vector<int>> data_max_idx(numberOfInputs_);

  bool do_min = false;
  bool do_sad = false;
  bool do_max = false;

  // Create diagrams for min, saddle and max persistence pairs
  for(int i = 0; i < numberOfInputs_; i++) {
    DiagramType &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      const auto &t = CTDiagram[j];

      const auto nt1 = std::get<1>(t);
      const auto nt2 = std::get<3>(t);
      const auto dt = std::get<4>(t);

      if(dt > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          data_max[i].push_back(t);
          data_max_idx[i].push_back(j);
          do_max = true;
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            data_max[i].push_back(t);
            data_max_idx[i].push_back(j);
            do_max = true;
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            data_min[i].push_back(t);
            data_min_idx[i].push_back(j);
            do_min = true;
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            data_sad[i].push_back(t);
            data_sad_idx[i].push_back(j);
            do_sad = true;
          }
        }
      }
    }
  }

  DiagramType barycenter_min{}, barycenter_sad{}, barycenter_max{};
  std::vector<std::vector<MatchingType>> matching_min{}, matching_sad{},
    matching_max{};

  double total_cost = 0;
  if(do_min && do_max) {
    time_limit_ = time_limit_ / 2;
  }
  if(do_sad) {
    time_limit_ = time_limit_ / 3;
  }

  const auto getRunner = [this]() -> PDBarycenter {
    PDBarycenter runner{};
    runner.setDebugLevel(this->debugLevel_);
    runner.setThreadNumber(this->threadNumber_);
    runner.setWasserstein(this->wasserstein_);
    runner.setNumberOfInputs(this->numberOfInputs_);
    runner.setUseProgressive(this->use_progressive_);
    runner.setTimeLimit(this->time_limit_);
    runner.setGeometricalFactor(this->alpha_);
    runner.setDeterministic(this->deterministic_);
    runner.setLambda(this->lambda_);
    runner.setMethod(this->method_);
    runner.setEarlyStoppage(this->early_stoppage_);
    runner.setEpsilonDecreases(this->epsilon_decreases_);
    runner.setReinitPrices(this->reinit_prices_);
    return runner;
  };

  if(do_min) {
    printMsg("Computing Minima barycenter...");
    PDBarycenter bary_min = getRunner();
    matching_min = bary_min.execute(data_min, barycenter_min);
    total_cost += bary_min.getCost();
  }

  if(do_sad) {
    printMsg("Computing Saddles barycenter...");
    PDBarycenter bary_sad = getRunner();
    matching_sad = bary_sad.execute(data_sad, barycenter_sad);
    total_cost += bary_sad.getCost();
  }

  if(do_max) {
    printMsg("Computing Maxima barycenter...");
    PDBarycenter bary_max = getRunner();
    matching_max = bary_max.execute(data_max, barycenter_max);
    total_cost += bary_max.getCost();
  }

  // Reconstruct matchings
  all_matchings.resize(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {

    if(do_min) {
      for(unsigned int j = 0; j < matching_min[i].size(); j++) {
        MatchingType t = matching_min[i][j];
        int bidder_id = std::get<0>(t);
        std::get<0>(t) = data_min_idx[i][bidder_id];
        if(std::get<1>(t) < 0) {
          std::get<1>(t) = -1;
        }
        all_matchings[i].push_back(t);
      }
    }

    if(do_sad) {
      for(unsigned int j = 0; j < matching_sad[i].size(); j++) {
        MatchingType t = matching_sad[i][j];
        int bidder_id = std::get<0>(t);
        std::get<0>(t) = data_sad_idx[i][bidder_id];
        if(std::get<1>(t) >= 0) {
          std::get<1>(t) = std::get<1>(t) + barycenter_min.size();
        } else {
          std::get<1>(t) = -1;
        }
        all_matchings[i].push_back(t);
      }
    }

    if(do_max) {
      for(unsigned int j = 0; j < matching_max[i].size(); j++) {
        MatchingType t = matching_max[i][j];
        int bidder_id = std::get<0>(t);
        std::get<0>(t) = data_max_idx[i][bidder_id];
        if(std::get<1>(t) >= 0) {
          std::get<1>(t)
            = std::get<1>(t) + barycenter_min.size() + barycenter_sad.size();
        } else {
          std::get<1>(t) = -1;
        }
        all_matchings[i].push_back(t);
      }
    }
  }
  // Reconstruct barycenter
  for(const auto &pair : barycenter_min) {
    barycenter.emplace_back(pair);
  }
  for(const auto &pair : barycenter_sad) {
    barycenter.emplace_back(pair);
  }
  for(const auto &pair : barycenter_max) {
    barycenter.emplace_back(pair);
  }

  // Recreate 3D critical coordinates of barycentric points
  std::vector<int> number_of_matchings_for_point(barycenter.size());
  std::vector<std::array<float, 3>> cords_1(barycenter.size());
  std::vector<std::array<float, 3>> cords_2(barycenter.size());

  for(size_t i = 0; i < all_matchings.size(); i++) {
    DiagramType &CTDiagram = intermediateDiagrams[i];
    for(size_t j = 0; j < all_matchings[i].size(); j++) {
      const auto &t = all_matchings[i][j];
      const int bidder_id = std::get<0>(t);
      const int bary_id = std::get<1>(t);

      const auto &bidder = CTDiagram[bidder_id];
      number_of_matchings_for_point[bary_id]++;
      cords_1[bary_id][0] += std::get<7>(bidder);
      cords_1[bary_id][1] += std::get<8>(bidder);
      cords_1[bary_id][2] += std::get<9>(bidder);
      cords_2[bary_id][0] += std::get<11>(bidder);
      cords_2[bary_id][1] += std::get<12>(bidder);
      cords_2[bary_id][2] += std::get<13>(bidder);
    }
  }

  for(size_t i = 0; i < barycenter.size(); i++) {
    if(number_of_matchings_for_point[i] > 0) {
      std::get<7>(barycenter[i])
        = cords_1[i][0] / number_of_matchings_for_point[i];
      std::get<8>(barycenter[i])
        = cords_1[i][1] / number_of_matchings_for_point[i];
      std::get<9>(barycenter[i])
        = cords_1[i][2] / number_of_matchings_for_point[i];
      std::get<11>(barycenter[i])
        = cords_2[i][0] / number_of_matchings_for_point[i];
      std::get<12>(barycenter[i])
        = cords_2[i][1] / number_of_matchings_for_point[i];
      std::get<13>(barycenter[i])
        = cords_2[i][2] / number_of_matchings_for_point[i];
    }
  }

  printMsg("Total cost : " + std::to_string(total_cost));
  printMsg("Complete", 1, tm.getElapsedTime(), threadNumber_);
}

using ttk::MatchingType;
using ttk::PDBarycenter;

std::vector<std::vector<MatchingType>>
  PDBarycenter::execute(const std::vector<DiagramType> &inputDiagrams,
                        ttk::DiagramType &barycenter) {
  return executeAuctionBarycenter(inputDiagrams, barycenter);
}

std::vector<std::vector<MatchingType>> PDBarycenter::executeAuctionBarycenter(
  const std::vector<DiagramType> &inputDiagrams, DiagramType &barycenter) {
  double total_time = 0;

  std::vector<std::vector<MatchingType>> previous_matchings;
  double min_persistence = 0;
  double min_cost = std::numeric_limits<double>::max();
  int last_min_cost_obtained = 0;

  const auto diagramType = std::get<5>(inputDiagrams[0][0]);
  const auto nt1 = std::get<1>(inputDiagrams[0][0]);
  const auto nt2 = std::get<3>(inputDiagrams[0][0]);

  this->setBidderDiagrams(inputDiagrams);
  this->setInitialBarycenter(inputDiagrams, min_persistence);

  double max_persistence = getMaxPersistence();

  std::vector<double> min_diag_price(numberOfInputs_);
  std::vector<double> min_price(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {
    min_diag_price[i] = 0;
    min_price[i] = 0;
  }

  int min_points_to_add = std::numeric_limits<int>::max();
  min_persistence = this->enrichCurrentBidderDiagrams(
    2 * max_persistence, min_persistence, min_diag_price, min_price,
    min_points_to_add, false);

  int n_iterations = 0;

  bool converged = false;
  bool finished = false;
  double total_cost;

  while(!finished) {
    Timer tm;

    n_iterations += 1;

    std::pair<std::unique_ptr<KDTree<double>>, std::vector<KDTree<double> *>>
      pair;
    bool use_kdt = false;
    // If the barycenter is empty, do not compute the kdt (or it will crash :/)
    // TODO Fix KDTree to handle empty inputs...
    if(barycenter_goods_[0].size() > 0) {
      pair = this->getKDTree();
      use_kdt = true;
    }

    std::vector<std::vector<MatchingType>> all_matchings(numberOfInputs_);
    std::vector<int> sizes(numberOfInputs_);
    for(int i = 0; i < numberOfInputs_; i++) {
      sizes[i] = current_bidder_diagrams_[i].size();
    }

    total_cost = 0;

    barycenter.clear();
    for(size_t j = 0; j < barycenter_goods_[0].size(); j++) {
      Good &g = barycenter_goods_[0][j];
      const auto t = std::make_tuple(0, nt1, 0, nt2, g.getPersistence(),
                                     diagramType, g.x_, 0, 0, 0, g.y_, 0, 0, 0);
      barycenter.push_back(t);
    }

    runMatchingAuction(total_cost, sizes, *pair.first, pair.second,
                       min_diag_price, all_matchings, use_kdt);

    this->printMsg("Barycenter cost : " + std::to_string(total_cost),
                   debug::Priority::DETAIL);

    if(converged) {
      finished = true;
    }

    if(!finished) {
      updateBarycenter(all_matchings);

      if(min_cost > total_cost) {
        min_cost = total_cost;
        last_min_cost_obtained = 0;
      } else {
        last_min_cost_obtained += 1;
      }

      converged = converged || last_min_cost_obtained > 1;
    }

    previous_matchings = std::move(all_matchings);
    // END OF TIMER
    total_time += tm.getElapsedTime();

    for(unsigned int i = 0; i < barycenter_goods_.size(); ++i) {
      for(size_t j = 0; j < barycenter_goods_[i].size(); ++j) {
        barycenter_goods_[i][j].setPrice(0);
      }
    }
    for(unsigned int i = 0; i < current_bidder_diagrams_.size(); ++i) {
      for(size_t j = 0; j < current_bidder_diagrams_[i].size(); ++j) {
        current_bidder_diagrams_[i][j].setDiagonalPrice(0);
      }
    }
    for(int i = 0; i < numberOfInputs_; i++) {
      min_diag_price[i] = 0;
      min_price[i] = 0;
    }
  }
  barycenter.clear();
  for(size_t j = 0; j < barycenter_goods_[0].size(); j++) {
    Good &g = barycenter_goods_[0][j];
    const auto t = std::make_tuple(0, nt1, 0, nt2, g.getPersistence(),
                                   diagramType, g.x_, 0, 0, 0, g.y_, 0, 0, 0);
    barycenter.push_back(t);
  }

  cost_ = sqrt(total_cost);
  std::vector<std::vector<MatchingType>> corrected_matchings
    = correctMatchings(previous_matchings);
  return corrected_matchings;
}

void PDBarycenter::runMatching(
  double &total_cost,
  const double epsilon,
  const std::vector<int> &sizes,
  KDTree<double> &kdt,
  std::vector<KDTree<double> *> &correspondance_kdt_map,
  std::vector<double> &min_diag_price,
  std::vector<double> &min_price,
  std::vector<std::vector<MatchingType>> &all_matchings,
  const bool use_kdt,
  const int actual_distance) {

  Timer time_matchings;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic, 1)
#endif
  for(int i = 0; i < numberOfInputs_; i++) {
    PersistenceDiagramAuction auction(
      current_bidder_diagrams_[i], barycenter_goods_[i], wasserstein_,
      geometrical_factor_, lambda_, 0.01, kdt, correspondance_kdt_map, epsilon,
      min_diag_price[i], use_kdt);
    int n_biddings = 0;
    auction.buildUnassignedBidders();
    auction.reinitializeGoods();
    auction.runAuctionRound(n_biddings, i);
    auction.updateDiagonalPrices();
    min_diag_price[i] = auction.getMinimalDiagonalPrice();
    min_price[i] = getMinimalPrice(i);
    std::vector<MatchingType> matchings;
    double cost = auction.getMatchingsAndDistance(matchings, true);
    all_matchings[i] = matchings;
    if(actual_distance) {
      total_cost += cost;
    } else {
      total_cost += cost * cost;
    }

    const double quotient
      = epsilon * auction.getAugmentedNumberOfBidders() / cost;
    precision_[i] = quotient < 1 ? 1. / sqrt(1 - quotient) - 1 : 10;
    current_bidder_diagrams_[i].resize(sizes[i]);
  }
}

void PDBarycenter::runMatchingAuction(
  double &total_cost,
  const std::vector<int> &sizes,
  KDTree<double> &kdt,
  std::vector<KDTree<double> *> &correspondance_kdt_map,
  std::vector<double> &min_diag_price,
  std::vector<std::vector<MatchingType>> &all_matchings,
  const bool use_kdt) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic, 1)
#endif
  for(int i = 0; i < numberOfInputs_; i++) {
    PersistenceDiagramAuction auction(
      current_bidder_diagrams_[i], barycenter_goods_[i], wasserstein_,
      geometrical_factor_, lambda_, 0.01, kdt, correspondance_kdt_map,
      min_diag_price[i], use_kdt);
    std::vector<MatchingType> matchings;
    double cost = auction.run(matchings);
    all_matchings[i] = matchings;

    total_cost += cost * cost;
    current_bidder_diagrams_[i].resize(sizes[i]);
  }
}

bool PDBarycenter::hasBarycenterConverged(
  std::vector<std::vector<MatchingType>> &matchings,
  std::vector<std::vector<MatchingType>> &previous_matchings) {

  if(points_added_ > 0 || points_deleted_ > 0
     || previous_matchings.size() == 0) {
    return false;
  }

  for(unsigned int j = 0; j < matchings.size(); j++) {
    for(unsigned int i = 0; i < matchings[j].size(); i++) {
      MatchingType t = matchings[j][i];
      MatchingType previous_t = previous_matchings[j][i];

      if(std::get<1>(t) != std::get<1>(previous_t)
         && (std::get<0>(t) >= 0 && std::get<0>(previous_t) >= 0)) {
        return false;
      }
    }
  }
  return true;
}

std::vector<std::vector<MatchingType>> PDBarycenter::correctMatchings(
  std::vector<std::vector<MatchingType>> previous_matchings) {

  std::vector<std::vector<MatchingType>> corrected_matchings(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {
    // 1. Invert the current_bidder_ids_ vector
    std::vector<int> new_to_old_id(current_bidder_diagrams_[i].size());
    for(unsigned int j = 0; j < current_bidder_ids_[i].size(); j++) {
      int new_id = current_bidder_ids_[i][j];
      if(new_id >= 0) {
        new_to_old_id[new_id] = j;
      }
    }
    // 2. Reconstruct the matchings
    std::vector<MatchingType> matchings_diagram_i;
    for(unsigned int j = 0; j < previous_matchings[i].size(); j++) {
      MatchingType m = previous_matchings[i][j];
      int new_id = std::get<0>(m);
      if(new_id >= 0 && std::get<1>(m) >= 0) {
        std::get<0>(m) = new_to_old_id[new_id];
        matchings_diagram_i.push_back(m);
      }
    }
    corrected_matchings[i] = matchings_diagram_i;
  }
  return corrected_matchings;
}

double PDBarycenter::updateBarycenter(
  std::vector<std::vector<MatchingType>> &matchings) {
  // 1. Initialize variables used in the sequel
  Timer t_update;
  unsigned int n_goods = barycenter_goods_[0].size();

  unsigned int n_diagrams = current_bidder_diagrams_.size();
  points_added_ = 0;
  points_deleted_ = 0;
  double max_shift = 0;

  std::vector<unsigned int> count_diag_matchings(
    n_goods); // Number of diagonal matchings for each point of the barycenter
  std::vector<double> x(n_goods);
  std::vector<double> y(n_goods);
  std::vector<double> crit_coords_x(n_goods);
  std::vector<double> crit_coords_y(n_goods);
  std::vector<double> crit_coords_z(n_goods);
  for(unsigned int i = 0; i < n_goods; i++) {
    count_diag_matchings[i] = 0;
    x[i] = 0;
    y[i] = 0;
    crit_coords_x[i] = 0;
    crit_coords_y[i] = 0;
    crit_coords_z[i] = 0;
  }
  std::vector<double> min_prices(n_diagrams);
  for(unsigned int j = 0; j < n_diagrams; j++) {
    min_prices[j] = std::numeric_limits<double>::max();
  }

  std::vector<Bidder *>
    points_to_append; // Will collect bidders linked to diagonal
  // 2. Preprocess the matchings
  for(unsigned int j = 0; j < matchings.size(); j++) {
    for(unsigned int i = 0; i < matchings[j].size(); i++) {
      int bidder_id = std::get<0>(matchings[j][i]);
      int good_id = std::get<1>(matchings[j][i]);
      if(good_id < 0 && bidder_id >= 0) {
        // Future new barycenter point
        points_to_append.push_back(&current_bidder_diagrams_[j][bidder_id]);
      }

      else if(good_id >= 0 && bidder_id >= 0) {
        // Update coordinates (to be divided by the number of diagrams later on)
        x[good_id] += current_bidder_diagrams_[j][bidder_id].x_;
        y[good_id] += current_bidder_diagrams_[j][bidder_id].y_;
        if(geometrical_factor_ < 1) {
          const auto critical_coordinates
            = current_bidder_diagrams_[j][bidder_id].GetCriticalCoordinates();
          crit_coords_x[good_id] += critical_coordinates[0];
          crit_coords_y[good_id] += critical_coordinates[1];
          crit_coords_z[good_id] += critical_coordinates[2];
        }
      } else if(good_id >= 0 && bidder_id < 0) {
        // Counting the number of times this barycenter point is linked to the
        // diagonal
        count_diag_matchings[good_id] = count_diag_matchings[good_id] + 1;
      }
    }
  }

  // 3. Update the previous points of the barycenter
  for(unsigned int i = 0; i < n_goods; i++) {
    if(count_diag_matchings[i] < n_diagrams) {
      // Barycenter point i is matched at least to one off-diagonal bidder
      // 3.1 Compute the arithmetic mean of off-diagonal bidders linked to it
      double x_bar = x[i] / (double)(n_diagrams - count_diag_matchings[i]);
      double y_bar = y[i] / (double)(n_diagrams - count_diag_matchings[i]);
      // 3.2 Compute the new coordinates of the point (the more linked to the
      // diagonal it was, the closer to the diagonal it'll be)
      double new_x = ((double)(n_diagrams - count_diag_matchings[i]) * x_bar
                      + (double)count_diag_matchings[i] * (x_bar + y_bar) / 2.)
                     / (double)n_diagrams;
      double new_y = ((double)(n_diagrams - count_diag_matchings[i]) * y_bar
                      + (double)count_diag_matchings[i] * (x_bar + y_bar) / 2.)
                     / (double)n_diagrams;
      // TODO Weight by persistence
      double new_crit_coord_x
        = crit_coords_x[i] / (double)(n_diagrams - count_diag_matchings[i]);
      double new_crit_coord_y
        = crit_coords_y[i] / (double)(n_diagrams - count_diag_matchings[i]);
      double new_crit_coord_z
        = crit_coords_z[i] / (double)(n_diagrams - count_diag_matchings[i]);

      // 3.3 Compute and store how much the point has shifted
      // TODO adjust shift with geometrical_factor_
      double dx = barycenter_goods_[0][i].x_ - new_x;
      double dy = barycenter_goods_[0][i].y_ - new_y;
      double shift = Geometry::pow(std::abs(dx), wasserstein_)
                     + Geometry::pow(std::abs(dy), wasserstein_);
      if(shift > max_shift) {
        max_shift = shift;
      }
      // 3.4 Update the position of the point
      for(unsigned int j = 0; j < n_diagrams; j++) {
        barycenter_goods_[j][i].SetCoordinates(new_x, new_y);
        if(geometrical_factor_ < 1) {
          barycenter_goods_[j][i].SetCriticalCoordinates(
            new_crit_coord_x, new_crit_coord_y, new_crit_coord_z);
        }
        if(barycenter_goods_[j][i].getPrice() < min_prices[j]) {
          min_prices[j] = barycenter_goods_[j][i].getPrice();
        }
      }
      // TODO Reinitialize/play with prices here if you wish
    }
  }
  for(unsigned int j = 0; j < n_diagrams; j++) {
    if(min_prices[j] >= std::numeric_limits<double>::max() / 2.) {
      min_prices[j] = 0;
    }
  }

  // 4. Delete off-diagonal barycenter points not linked to any
  // off-diagonal bidder
  for(unsigned int i = 0; i < n_goods; i++) {
    if(count_diag_matchings[i] == n_diagrams) {
      points_deleted_ += 1;
      double shift
        = 2
          * Geometry::pow(
            barycenter_goods_[0][i].getPersistence() / 2., wasserstein_);
      if(shift > max_shift) {
        max_shift = shift;
      }
      for(unsigned int j = 0; j < n_diagrams; j++) {
        barycenter_goods_[j][i].id_ = -1;
      }
    }
  }

  // 5. Append the new points to the barycenter
  for(unsigned int k = 0; k < points_to_append.size(); k++) {
    points_added_ += 1;
    Bidder *b = points_to_append[k];
    double gx
      = (b->x_ + (n_diagrams - 1) * (b->x_ + b->y_) / 2.) / (n_diagrams);
    double gy
      = (b->y_ + (n_diagrams - 1) * (b->x_ + b->y_) / 2.) / (n_diagrams);
    const auto critical_coordinates = b->GetCriticalCoordinates();
    for(unsigned int j = 0; j < n_diagrams; j++) {
      Good g(gx, gy, false, barycenter_goods_[j].size());
      g.setPrice(min_prices[j]);
      if(geometrical_factor_ < 1) {
        g.SetCriticalCoordinates(critical_coordinates);
      }
      barycenter_goods_[j].emplace_back(g);
      double shift
        = 2
          * Geometry::pow(
            barycenter_goods_[j][g.id_].getPersistence() / 2., wasserstein_);
      if(shift > max_shift) {
        max_shift = shift;
      }
    }
  }

  // 6. Finally, recreate barycenter_goods
  for(unsigned int j = 0; j < n_diagrams; j++) {
    int count = 0;
    GoodDiagram new_barycenter{};
    for(size_t i = 0; i < barycenter_goods_[j].size(); i++) {
      Good g = barycenter_goods_[j][i];
      if(g.id_ != -1) {
        g.id_ = count;
        new_barycenter.emplace_back(g);
        count++;
      }
    }
    barycenter_goods_[j] = new_barycenter;
  }

  return max_shift;
}

void PDBarycenter::setBidderDiagrams(
  const std::vector<DiagramType> &inputDiagrams) {

  for(int i = 0; i < numberOfInputs_; i++) {
    const auto &CTDiagram = inputDiagrams[i];

    BidderDiagram bidders;
    for(unsigned int j = 0; j < CTDiagram.size(); j++) {
      // Add bidder to bidders
      Bidder b(CTDiagram[j], j, lambda_);

      b.setPositionInAuction(bidders.size());
      bidders.emplace_back(b);
      if(b.isDiagonal() || b.x_ == b.y_) {
        this->printWrn("Diagonal point in diagram !!!");
      }
    }
    bidder_diagrams_.push_back(bidders);
    current_bidder_diagrams_.emplace_back(BidderDiagram{});
    std::vector<int> ids(bidders.size());
    for(unsigned int j = 0; j < ids.size(); j++) {
      ids[j] = -1;
    }
    current_bidder_ids_.push_back(ids);
  }
  return;
}

double PDBarycenter::enrichCurrentBidderDiagrams(
  double previous_min_persistence,
  double min_persistence,
  std::vector<double> initial_diagonal_prices,
  std::vector<double> initial_off_diagonal_prices,
  int min_points_to_add,
  bool add_points_to_barycenter) {

  double new_min_persistence = min_persistence;

  // 1. Get size of the largest current diagram, deduce the maximal number of
  // points to append
  size_t max_diagram_size = 0;
  for(int i = 0; i < numberOfInputs_; i++) {
    if(current_bidder_diagrams_[i].size() > max_diagram_size) {
      max_diagram_size = current_bidder_diagrams_[i].size();
    }
  }
  int max_points_to_add = std::max(
    min_points_to_add, min_points_to_add + (int)(max_diagram_size / 10));

  // 2. Get which points can be added, deduce the new minimal persistence
  std::vector<std::vector<int>> candidates_to_be_added(numberOfInputs_);
  std::vector<std::vector<int>> idx(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {

    std::vector<double> persistences;
    for(size_t j = 0; j < bidder_diagrams_[i].size(); j++) {
      Bidder b = bidder_diagrams_[i][j];
      double persistence = b.getPersistence();
      if(persistence >= min_persistence
         && persistence < previous_min_persistence) {
        candidates_to_be_added[i].push_back(j);
        idx[i].push_back(idx[i].size());
        persistences.push_back(persistence);
      }
    }
    sort(idx[i].begin(), idx[i].end(), [&persistences](int &a, int &b) {
      return ((persistences[a] > persistences[b])
              || ((persistences[a] == persistences[b]) && (a > b)));
    });
    int size = candidates_to_be_added[i].size();
    if(size >= max_points_to_add) {
      double last_persistence_added
        = persistences[idx[i][max_points_to_add - 1]];
      if(last_persistence_added > new_min_persistence) {
        new_min_persistence = last_persistence_added;
      }
    }
  }

  // 3. Add the points to the current diagrams

  // only to give determinism
  int compteur_for_adding_points = 0;

  for(int i = 0; i < numberOfInputs_; i++) {
    int size = candidates_to_be_added[i].size();
    for(int j = 0; j < std::min(max_points_to_add, size); j++) {
      Bidder b = bidder_diagrams_[i][candidates_to_be_added[i][idx[i][j]]];
      if(b.getPersistence() >= new_min_persistence) {
        b.id_ = current_bidder_diagrams_[i].size();
        b.setPositionInAuction(current_bidder_diagrams_[i].size());
        b.setDiagonalPrice(initial_diagonal_prices[i]);
        current_bidder_diagrams_[i].emplace_back(b);
        // b.id_ --> position of b in current_bidder_diagrams_[i]
        current_bidder_ids_[i][candidates_to_be_added[i][idx[i][j]]]
          = current_bidder_diagrams_[i].size() - 1;

        int to_be_added_to_barycenter
          = deterministic_ ? compteur_for_adding_points % numberOfInputs_
                           : rand() % numberOfInputs_;
        // We add the bidder as a good with probability 1/n_diagrams
        if(to_be_added_to_barycenter == 0 && add_points_to_barycenter) {
          for(int k = 0; k < numberOfInputs_; k++) {
            Good g(b.x_, b.y_, false, barycenter_goods_[k].size());
            g.setPrice(initial_off_diagonal_prices[k]);
            g.SetCriticalCoordinates(b.coords_);
            barycenter_goods_[k].emplace_back(g);
          }
        }
      }
      compteur_for_adding_points++;
    }
  }
  return new_min_persistence;
}

double PDBarycenter::getMaxPersistence() {
  double max_persistence = 0;
  for(int i = 0; i < numberOfInputs_; i++) {
    BidderDiagram &D = bidder_diagrams_[i];
    for(size_t j = 0; j < D.size(); j++) {
      // Add bidder to bidders
      Bidder &b = D[j];
      double persistence = b.getPersistence();
      if(persistence > max_persistence) {
        max_persistence = persistence;
      }
    }
  }
  return max_persistence;
}

double PDBarycenter::getMinimalPrice(int i) {
  double min_price = std::numeric_limits<double>::max();

  GoodDiagram &D = barycenter_goods_[i];
  if(D.size() == 0) {
    return 0;
  }
  for(size_t j = 0; j < D.size(); j++) {
    Good &b = D[j];
    double price = b.getPrice();
    if(price < min_price) {
      min_price = price;
    }
  }
  if(min_price >= std::numeric_limits<double>::max() / 2.) {
    return 0;
  }
  return min_price;
}

double PDBarycenter::getLowestPersistence() {
  double lowest_persistence = std::numeric_limits<double>::max();
  for(int i = 0; i < numberOfInputs_; i++) {
    BidderDiagram &D = bidder_diagrams_[i];
    for(size_t j = 0; j < D.size(); j++) {
      // Add bidder to bidders
      Bidder &b = D[j];
      double persistence = b.getPersistence();
      if(persistence < lowest_persistence && persistence > 0) {
        lowest_persistence = persistence;
      }
    }
  }
  if(lowest_persistence >= std::numeric_limits<double>::max() / 2.) {
    return 0;
  }
  return lowest_persistence;
}

void PDBarycenter::setInitialBarycenter(
  const std::vector<DiagramType> &inputDiagrams, double min_persistence) {
  int size = 0;
  int random_idx;
  int iter = 0;
  while(size == 0) {
    random_idx
      = deterministic_ ? iter % numberOfInputs_ : rand() % numberOfInputs_;
    const auto &CTDiagram = inputDiagrams[random_idx];
    size = CTDiagram.size();
    for(int i = 0; i < numberOfInputs_; i++) {
      GoodDiagram goods;
      int count = 0;
      for(unsigned int j = 0; j < CTDiagram.size(); j++) {
        // Add bidder to bidders
        Good g(CTDiagram[j], count, lambda_);
        if(g.getPersistence() >= min_persistence) {
          goods.emplace_back(g);
          count++;
        }
      }
      if(barycenter_goods_.size() < (unsigned int)(i + 1)) {
        barycenter_goods_.push_back(goods);
      } else {
        barycenter_goods_[i] = goods;
      }
    }
    size = barycenter_goods_[0].size();
    iter++;
  }
}

typename PDBarycenter::KDTreePair PDBarycenter::getKDTree() const {
  Timer tm;
  auto kdt
    = std::unique_ptr<KDTree<double>>(new KDTree<double>{true, wasserstein_});

  const int dimension = geometrical_factor_ >= 1 ? 2 : 5;

  std::vector<double> coordinates;
  std::vector<std::vector<double>> weights;

  for(size_t i = 0; i < barycenter_goods_[0].size(); i++) {
    const Good &g = barycenter_goods_[0][i];
    coordinates.push_back(geometrical_factor_ * g.x_);
    coordinates.push_back(geometrical_factor_ * g.y_);
    if(geometrical_factor_ < 1) {
      coordinates.push_back((1 - geometrical_factor_) * g.coords_[0]);
      coordinates.push_back((1 - geometrical_factor_) * g.coords_[1]);
      coordinates.push_back((1 - geometrical_factor_) * g.coords_[2]);
    }
  }

  for(unsigned int idx = 0; idx < barycenter_goods_.size(); idx++) {
    std::vector<double> empty_weights;
    weights.push_back(empty_weights);
    for(size_t i = 0; i < barycenter_goods_[idx].size(); i++) {
      const Good &g = barycenter_goods_[idx][i];
      weights[idx].push_back(g.getPrice());
    }
  }
  // Correspondance map : position in barycenter_goods_ --> KDT node

  auto correspondance_kdt_map
    = kdt->build(coordinates.data(), barycenter_goods_[0].size(), dimension,
                 weights, barycenter_goods_.size());
  this->printMsg(" Building KDTree", 1, tm.getElapsedTime(),
                 debug::LineMode::NEW, debug::Priority::VERBOSE);
  return std::make_pair(std::move(kdt), correspondance_kdt_map);
}

double PDBarycenter::computeRealCost() {
  double total_real_cost = 0;
  for(int i = 0; i < numberOfInputs_; i++) {
    PersistenceDiagramAuction auction(
      wasserstein_, geometrical_factor_, lambda_, 0.01, true);
    GoodDiagram current_barycenter = barycenter_goods_[0];
    BidderDiagram current_bidder_diagram = bidder_diagrams_[i];
    auction.BuildAuctionDiagrams(current_bidder_diagram, current_barycenter);
    double cost = auction.run();
    total_real_cost += cost * cost;
  }
  return sqrt(total_real_cost);
}

bool PDBarycenter::isPrecisionObjectiveMet(double precision_objective,
                                           int mode) {
  if(mode == 0) { // ABSOLUTE PRECISION
    for(int i_input = 0; i_input < numberOfInputs_; i_input++) {
      if(precision_[i_input] > precision_objective) {
        return false;
      }
    }
  } else if(mode == 1) { // AVERAGE PRECISION
    double average_precision
      = std::accumulate(precision_.begin(), precision_.end(), 0.0)
        / numberOfInputs_;
    if(average_precision > precision_objective) {
      return false;
    }
  }
  return true;
}
