/// \ingroup base
/// \class ttk::PDBarycenter
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering

#pragma once

#include "DataTypes.h"
#include <KDTree.h>
#include <PersistenceDiagramAuction.h>
#include <PersistenceDiagramBarycenter.h>

#include <limits>

namespace ttk {

  class PDBarycenter : public Debug {
    enum class ComputeMethod { PARTIAL_BIDDING, MUNKRES, AUCTION };

  public:
    PDBarycenter() {
      this->setDebugMsgPrefix("PersistenceDiagramBarycenter");
    }

    std::vector<std::vector<MatchingType>> execute(DiagramType &barycenter);
    std::vector<std::vector<MatchingType>>
      executeMunkresBarycenter(DiagramType &barycenter);
    std::vector<std::vector<MatchingType>>
      executeAuctionBarycenter(DiagramType &barycenter);
    std::vector<std::vector<MatchingType>>
      executePartialBiddingBarycenter(DiagramType &barycenter);

    void setBidderDiagrams();
    double
      enrichCurrentBidderDiagrams(double previous_min_persistence,
                                  double min_persistence,
                                  std::vector<double> initial_diagonal_prices,
                                  std::vector<double> initial_prices,
                                  int min_points_to_add,
                                  bool add_points_to_barycenter = true);
    void setInitialBarycenter(double min_persistence);
    double getMaxPersistence();
    double getLowestPersistence();
    double getMinimalPrice(int i);
    using KDTreePair = std::pair<typename KDTree<double>::KDTreeRoot,
                                 typename KDTree<double>::KDTreeMap>;
    KDTreePair getKDTree() const;

    void runMatching(double *total_cost,
                     double epsilon,
                     std::vector<int> sizes,
                     KDTree<double> &kdt,
                     std::vector<KDTree<double> *> &correspondance_kdt_map,
                     std::vector<double> *min_diag_price,
                     std::vector<double> *min_price,
                     std::vector<std::vector<MatchingType>> *all_matchings,
                     bool use_kdt,
                     int compute_only_distance);

    void
      runMatchingAuction(double *total_cost,
                         std::vector<int> sizes,
                         KDTree<double> &kdt,
                         std::vector<KDTree<double> *> &correspondance_kdt_map,
                         std::vector<double> *min_diag_price,
                         std::vector<std::vector<MatchingType>> *all_matchings,
                         bool use_kdt);

    double updateBarycenter(std::vector<std::vector<MatchingType>> &matchings);

    double computeRealCost();
    bool isPrecisionObjectiveMet(double, int);
    bool hasBarycenterConverged(
      std::vector<std::vector<MatchingType>> &matchings,
      std::vector<std::vector<MatchingType>> &previous_matchings);

    std::vector<std::vector<MatchingType>> correctMatchings(
      std::vector<std::vector<MatchingType>> previous_matchings);

    bool is_matching_stable();

    inline double getEpsilon(double rho) const {
      return rho * rho / 8.0;
    }
    inline double getRho(double epsilon) const {
      return std::sqrt(8.0 * epsilon);
    }

    inline void setDeterministic(const bool deterministic) {
      deterministic_ = deterministic;
    }

    inline void setMethod(const int method) {
      if(method == 0) {
        method_ = ComputeMethod::PARTIAL_BIDDING;
      } else if(method == 1) {
        method_ = ComputeMethod::MUNKRES;
      } else if(method == 2) {
        method_ = ComputeMethod::AUCTION;
      }
    }

    inline void setDiagrams(std::vector<DiagramType> *const data) {
      inputDiagrams_ = data;
    }

    inline void setNumberOfInputs(const int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
      precision_.resize(numberOfInputs_);
    }

    inline void setWasserstein(const int wasserstein) {
      wasserstein_ = wasserstein;
    }

    inline void setUseProgressive(const bool use_progressive) {
      use_progressive_ = use_progressive;
    }

    inline void setTimeLimit(const double time_limit) {
      time_limit_ = time_limit;
    }

    inline void setGeometricalFactor(const double geometrical_factor) {
      geometrical_factor_ = geometrical_factor;
    }

    inline void setLambda(const double lambda) {
      lambda_ = lambda;
    }

    inline void setCurrentBidders(const std::vector<BidderDiagram> &diagrams) {
      current_bidder_diagrams_ = diagrams;
    }

    inline void
      setCurrentBarycenter(const std::vector<GoodDiagram> &barycenters) {
      barycenter_goods_ = barycenters;
    }

    inline std::vector<BidderDiagram> &getCurrentBidders() {
      return current_bidder_diagrams_;
    }

    inline std::vector<GoodDiagram> &getCurrentBarycenter() {
      return barycenter_goods_;
    }

    inline void setReinitPrices(const bool reinit_prices) {
      reinit_prices_ = reinit_prices;
    }

    inline void setEpsilonDecreases(const bool epsilon_decreases) {
      epsilon_decreases_ = epsilon_decreases;
    }

    inline void setEarlyStoppage(const bool early_stoppage) {
      early_stoppage_ = early_stoppage;
    }

    inline void setDiagramType(const int diagramType) {
      diagramType_ = diagramType;
      if(diagramType_ == 0) {
        nt1_ = CriticalType::Local_minimum;
        nt2_ = CriticalType::Saddle1;
      } else if(diagramType_ == 1) {
        nt1_ = CriticalType::Saddle1;
        nt2_ = CriticalType::Saddle2;
      } else {
        nt1_ = CriticalType::Saddle2;
        nt2_ = CriticalType::Local_maximum;
      }
    }

    double getCost() {
      return cost_;
    }

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

  protected:
    std::vector<double> precision_{};
    ComputeMethod method_{ComputeMethod::PARTIAL_BIDDING};
    int wasserstein_{2};

    double geometrical_factor_{1.0};
    // lambda_ : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda_ = 1 : extremum (min if pair min-sad, max
    // if pair sad-max) lambda_ = 0 : saddle (awful stability) lambda_ = 1/2 :
    // middle of the 2 critical points of the pair (bad stability)
    double lambda_{};

    int diagramType_{};
    CriticalType nt1_{};
    CriticalType nt2_{};
    double cost_{};
    int numberOfInputs_{};
    double time_limit_{std::numeric_limits<double>::max()};
    double epsilon_min_{1e-5};
    std::vector<DiagramType> *inputDiagrams_{};

    int points_added_{};
    int points_deleted_{};

    std::vector<std::vector<double>> all_matchings_{};
    std::vector<std::vector<double>> all_old_matchings_{};
    std::vector<BidderDiagram> bidder_diagrams_{};
    std::vector<BidderDiagram> current_bidder_diagrams_{};
    std::vector<std::vector<int>> current_bidder_ids_{};
    std::vector<GoodDiagram> barycenter_goods_{};

    bool deterministic_{false}; // to kill any randomness
    bool use_progressive_{true};
    bool reinit_prices_{true};
    bool epsilon_decreases_{true};
    bool early_stoppage_{true};
  };
} // namespace ttk
