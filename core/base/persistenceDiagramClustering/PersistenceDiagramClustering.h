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

#include <PDClustering.h>
#include <PersistenceDiagramBarycenter.h>

namespace ttk {

  class PersistenceDiagramClustering : virtual public Debug {

  public:
    PersistenceDiagramClustering() {
      this->setDebugMsgPrefix("PersistenceDiagramClustering");
    }

    std::vector<int> execute(
      std::vector<DiagramType> &intermediateDiagrams,
      std::vector<DiagramType> &centroids,
      std::vector<std::vector<std::vector<MatchingType>>> &all_matchings);

  protected:
    // Critical pairs used for clustering
    // 0:min-saddles ; 1:saddles-saddles ; 2:sad-max ; else : all

    int DistanceWritingOptions{0};
    int PairTypeClustering{-1};
    bool Deterministic{true};
    int WassersteinMetric{2};

    bool UseProgressive{true};

    bool ForceUseOfAlgorithm{false};
    bool UseInterruptible{true};
    double Alpha{1.0};
    bool UseAdditionalPrecision{false};
    double DeltaLim{0.01};
    double Lambda{1.0};
    double TimeLimit{999999};

    int NumberOfClusters{1};
    bool UseAccelerated{false};
    bool UseKmeansppInit{false};

    int points_added_;
    int points_deleted_;
  };

} // namespace ttk
