/// \ingroup base
/// \class ttk::PersistenceDiagramClustering
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \brief TTK processing package for the computation of Wasserstein barycenters
/// and K-Means clusterings of a set of persistence diagrams.
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa ttkPersistenceDiagramClustering

#pragma once

#include <PersistenceDiagramBarycenter.h>

namespace ttk {

  class PersistenceDiagramClustering : virtual public Debug {

  public:
    PersistenceDiagramClustering() {
      this->setDebugMsgPrefix("PersistenceDiagramClustering");
    }

    std::vector<int> executeClustering(
      std::vector<DiagramType> &intermediateDiagrams,
      std::vector<DiagramType> &centroids,
      std::vector<std::vector<std::vector<MatchingType>>> &all_matchings);

  protected:
    std::vector<int>
      run(std::vector<DiagramType> &final_centroids,
          std::vector<std::array<std::vector<std::vector<MatchingType>>, 3>>
            &all_matchings);

    double getMostPersistent(const int type = -1) const;
    double getLessPersistent(const int type = -1) const;

    inline const std::vector<std::vector<int>> &get_centroids_sizes() const {
      return this->centroids_sizes_;
    }

    std::vector<std::vector<double>> getMinDiagonalPrices() const;
    std::vector<std::vector<double>> getMinPrices() const;

    void correctMatchings(
      std::vector<std::array<std::vector<std::vector<MatchingType>>, 3>>
        &previous_matchings);

    inline double computeDistance(const BidderDiagram &D1,
                                  const BidderDiagram &D2,
                                  const double delta_lim) const {
      const auto D2_bis = diagramToCentroid(D2);
      return computeDistance(D1, D2_bis, delta_lim);
    }
    double computeDistance(const BidderDiagram &D1,
                           const GoodDiagram &D2,
                           const double delta_lim) const;
    double computeDistance(BidderDiagram *const D1,
                           const GoodDiagram *const D2,
                           const double delta_lim) const;
    inline double computeDistance(const GoodDiagram &D1,
                                  const GoodDiagram &D2,
                                  const double delta_lim) const {
      const auto D1_bis = centroidToDiagram(D1);
      return computeDistance(D1_bis, D2, delta_lim);
    }

    GoodDiagram centroidWithZeroPrices(const GoodDiagram &centroid) const;
    BidderDiagram centroidToDiagram(const GoodDiagram &centroid) const;
    GoodDiagram diagramToCentroid(const BidderDiagram &diagram) const;
    BidderDiagram diagramWithZeroPrices(const BidderDiagram &diagram) const;

    void setBidderDiagrams();
    inline void initializeEmptyClusters() {
      this->clustering_.clear();
      this->clustering_.resize(this->NumberOfClusters);
    }
    void initializeCentroids();
    void initializeCentroidsKMeanspp();
    void initializeAcceleratedKMeans();
    void initializeBarycenterComputers();
    void printDistancesToFile() const;
    void printMatchings(
      const std::vector<std::vector<std::vector<MatchingType>>> &) const;
    void printRealDistancesToFile() const;
    void printPricesToFile(const int) const;
    double computeRealCost() const;

    std::vector<double> enrichCurrentBidderDiagrams(
      const std::vector<double> &previous_min_persistence,
      const std::vector<double> &min_persistence,
      const std::vector<std::vector<double>> &initial_diagonal_prices,
      const std::vector<std::vector<double>> &initial_off_diagonal_points,
      const std::vector<int> &min_points_to_add,
      const bool add_points_to_barycenter,
      const bool first_enrichment);

    std::vector<std::vector<double>> getDistanceMatrix() const;
    void getCentroidDistanceMatrix();
    void computeDistanceToCentroid();

    void updateClusters();
    void invertClusters();
    void invertInverseClusters();
    void computeBarycenterForTwo(
      std::vector<std::array<std::vector<std::vector<MatchingType>>, 3>> &);

    void acceleratedUpdateClusters();
    std::vector<double> updateCentroidsPosition(
      std::vector<std::vector<double>> &min_price,
      std::vector<std::vector<double>> &min_diag_price,
      std::vector<std::array<std::vector<std::vector<MatchingType>>, 3>>
        &all_matchings,
      const int only_matchings);

    inline void resetDosToOriginalValues() {
      do_min_ = original_dos[0];
      do_sad_ = original_dos[1];
      do_max_ = original_dos[2];
    }

    inline void setDos(bool doMin, bool doSad, bool doMax) {
      do_min_ = doMin;
      do_sad_ = doSad;
      do_max_ = doMax;

      this->original_dos = {do_min_, do_sad_, do_max_};
    }

    inline void setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
    }
    inline void setNumberOfClusters(const int k) {
      NumberOfClusters = k;
    }
    inline void setWasserstein(const int &wasserstein) {
      WassersteinMetric = wasserstein;
    }
    inline void setUseProgressive(const bool use_progressive) {
      UseProgressive = use_progressive;
    }
    inline void setKMeanspp(const bool use_kmeanspp) {
      UseKmeansppInit = use_kmeanspp;
    }
    inline void setUseKDTree(const bool use_kdtree) {
      use_kdtree_ = use_kdtree;
    }
    inline void setAccelerated(const bool use_accelerated) {
      UseAccelerated = use_accelerated;
    }
    inline void setTimeLimit(const double time_limit) {
      TimeLimit = time_limit;
    }
    inline void setGeometricalFactor(const double geometrical_factor) {
      Alpha = geometrical_factor;
    }
    inline void setLambda(const double lambda) {
      Lambda = lambda;
    }
    inline void setForceUseOfAlgorithm(const bool forceUseOfAlgorithm) {
      ForceUseOfAlgorithm = forceUseOfAlgorithm;
    }
    inline void setDeterministic(const bool deterministic) {
      Deterministic = deterministic;
    }
    inline void setUseDeltaLim(const bool UseDeltaLim) {
      UseAdditionalPrecision = UseDeltaLim;
      epsilon_min_ = UseAdditionalPrecision ? 1e-8 : 5e-5;
    }
    inline void setDistanceWritingOptions(const int distanceWritingOptions) {
      DistanceWritingOptions = distanceWritingOptions;
    }
    inline void setDeltaLim(const double deltaLim) {
      DeltaLim = deltaLim;
    }

    inline void printClustering() const {
      std::string msg = "";
      for(int c = 0; c < NumberOfClusters; ++c) {
        msg.append(" Cluster " + std::to_string(c) + " = {");
        for(unsigned int idx = 0; idx < clustering_[c].size(); ++idx) {
          if(idx == clustering_[c].size() - 1) {
            msg.append(std::to_string(clustering_[c][idx]) + "}");
            this->printMsg(msg);
            msg = "";
          } else {
            msg.append(std::to_string(clustering_[c][idx]) + ", ");
          }
        }
      }
    }

    // data members
    std::vector<PDBarycenter> barycenter_computer_min_{};
    std::vector<PDBarycenter> barycenter_computer_sad_{};
    std::vector<PDBarycenter> barycenter_computer_max_{};

    bool barycenter_inputs_reset_flag{};
    bool precision_criterion_{false};
    bool precision_max_{false};
    bool precision_min_{false};
    bool precision_sad_{false};
    bool ForceUseOfAlgorithm{false};
    bool Deterministic{true};
    int WassersteinMetric{2};
    double Alpha{1.0};
    double DeltaLim{0.01};
    bool UseAdditionalPrecision{false};
    int DistanceWritingOptions{0};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double Lambda{1.0};

    int numberOfInputs_{};
    int NumberOfClusters{1};
    int PairTypeClustering{-1};
    bool UseProgressive{true};
    bool UseAccelerated{false};
    bool UseKmeansppInit{false};
    bool use_kdtree_{true};
    double TimeLimit{std::numeric_limits<double>::max()};

    double epsilon_min_{1e-8};
    std::array<double, 3> epsilon_{};
    double cost_{};
    double cost_min_{};
    double cost_sad_{};
    double cost_max_{};

    std::vector<std::vector<int>> current_bidder_ids_min_;
    std::vector<std::vector<int>> current_bidder_ids_sad_;
    std::vector<std::vector<int>> current_bidder_ids_max_;
    std::vector<DiagramType> inputDiagramsMin_;
    std::vector<DiagramType> inputDiagramsSaddle_;
    std::vector<DiagramType> inputDiagramsMax_;

    std::array<bool, 3> original_dos{false, false, false};

    bool do_min_{false};
    std::vector<BidderDiagram> bidder_diagrams_min_;
    std::vector<BidderDiagram> current_bidder_diagrams_min_;
    std::vector<GoodDiagram> centroids_min_;
    std::vector<GoodDiagram> centroids_with_price_min_;

    bool do_sad_{false};
    std::vector<BidderDiagram> bidder_diagrams_saddle_;
    std::vector<BidderDiagram> current_bidder_diagrams_saddle_;
    std::vector<GoodDiagram> centroids_saddle_;
    std::vector<GoodDiagram> centroids_with_price_saddle_;

    bool do_max_{false};
    std::vector<BidderDiagram> bidder_diagrams_max_;
    std::vector<BidderDiagram> current_bidder_diagrams_max_;
    std::vector<GoodDiagram> centroids_max_;
    std::vector<GoodDiagram> centroids_with_price_max_;

    std::vector<std::vector<int>> clustering_;
    std::vector<std::vector<int>> old_clustering_;
    std::vector<int> inv_clustering_;

    std::vector<std::vector<int>> centroids_sizes_;

    std::vector<bool> r_;
    std::vector<double> u_;
    std::vector<std::vector<double>> l_;
    std::vector<std::vector<double>> centroidsDistanceMatrix_{};
    std::vector<double> distanceToCentroid_{};

    int n_iterations_;
  };

} // namespace ttk
