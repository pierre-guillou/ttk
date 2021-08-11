/// \ingroup base
/// \class ttk::PersistenceDiagramBarycenter
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

#ifndef diagramTuple
#define diagramTuple                                                       \
  std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,            \
             ttk::CriticalType, dataType, ttk::SimplexId, dataType, float, \
             float, float, dataType, float, float, float>
#endif

#ifndef BNodeType
#define BNodeType ttk::CriticalType
#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2
#define BIdVertex ttk::SimplexId
#endif

// base code includes
#include <KDTree.h>
#include <PDBarycenter.h>
#include <PersistenceDiagramAuction.h>
#include <Wrapper.h>

#include <limits>

using namespace std;

namespace ttk {

  class PersistenceDiagramBarycenter : public Debug {
  public:
    PersistenceDiagramBarycenter() {
      this->setDebugMsgPrefix("PersistenceDiagramBarycenter");
    }

    void execute(
      std::vector<DiagramType> &intermediateDiagrams,
      DiagramType &barycenter,
      std::vector<std::vector<std::vector<MatchingType>>> &all_matchings);

    inline void setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
    }

    inline void setDeterministic(const bool deterministic) {
      deterministic_ = deterministic;
    }

    inline void setWasserstein(const std::string &wasserstein) {
      wasserstein_ = (wasserstein == "inf") ? -1 : stoi(wasserstein);
    }

    inline void setUseProgressive(const bool use_progressive) {
      if(use_progressive)
        epsilon_decreases_ = true;
      use_progressive_ = use_progressive;
    }

    inline void setAlpha(const double alpha) {
      alpha_ = alpha;
    }

    inline void setLambda(const double lambda) {
      lambda_ = lambda;
    }

    inline void setTimeLimit(const double time_limit) {
      time_limit_ = time_limit;
    }

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

    inline void setMethod(const int &method) {
      method_ = method;
    }

    inline void setReinitPrices(const bool reinit_prices) {
      reinit_prices_ = reinit_prices;
    }

    inline void setEpsilonDecreases(const bool epsilon_decreases) {
      if(use_progressive_)
        epsilon_decreases_ = true;
      else
        epsilon_decreases_ = epsilon_decreases;
    }

    inline void setEarlyStoppage(const bool early_stoppage) {
      early_stoppage_ = early_stoppage;
    }

  protected:
    int method_{};
    int wasserstein_{2};
    int numberOfInputs_{0};
    double alpha_{1.0};
    double lambda_{1.0};
    double time_limit_{1.0};

    int points_added_{};
    int points_deleted_{};

    std::vector<std::vector<double>> all_matchings_{};
    std::vector<std::vector<double>> all_old_matchings_{};
    std::vector<BidderDiagram> bidder_diagrams_{};
    std::vector<GoodDiagram> barycenter_goods_{};

    bool deterministic_{true};
    bool use_progressive_{true};
    bool reinit_prices_{true};
    bool epsilon_decreases_{true};
    bool early_stoppage_{};
  };

} // namespace ttk
