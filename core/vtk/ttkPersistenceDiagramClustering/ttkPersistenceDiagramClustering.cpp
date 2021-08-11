#include <ttkMacros.h>
#include <ttkPersistenceDiagramClustering.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>

using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramClustering);

ttkPersistenceDiagramClustering::ttkPersistenceDiagramClustering() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

int ttkPersistenceDiagramClustering::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else {
    return 0;
  }
  return 1;
}

int ttkPersistenceDiagramClustering::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else {
    return 0;
  }
  return 1;
}

void ttkPersistenceDiagramClustering::Modified() {
  needUpdate_ = true;
  ttkAlgorithm::Modified();
}

int ttkPersistenceDiagramClustering::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  Memory m;

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  // Flat storage for diagrams extracted from blocks
  std::vector<vtkUnstructuredGrid *> input;

  // Number of input diagrams
  int numInputs = 0;

  if(blocks != nullptr) {
    numInputs = blocks->GetNumberOfBlocks();
    input.resize(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      input[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
      if(this->GetMTime() < input[i]->GetMTime()) {
        needUpdate_ = true;
      }
    }
  }

  if(numInputs == 0) {
    this->printErr("No input detected");
    return 0;
  }

  // Get output pointers
  auto output_clusters = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_centroids = vtkMultiBlockDataSet::GetData(outputVector, 1);
  auto output_matchings = vtkMultiBlockDataSet::GetData(outputVector, 2);

  if(needUpdate_) {
    // clear data before computation
    intermediateDiagrams_ = {};
    all_matchings_ = {};
    final_centroids_ = {};

    intermediateDiagrams_.resize(numInputs);
    all_matchings_.resize(3);

    // store the persistence of every min-max global pair
    std::vector<double> max_persistences(numInputs);

    for(int i = 0; i < numInputs; i++) {
      auto &diag{this->intermediateDiagrams_[i]};
      const auto ret = VTUToDiagram(diag, input[i], *this);
      if(ret < 0) {
        this->printErr("Could not read Persistence Diagram");
        return 0;
      }
      if(this->NumberOfClusters > 1) {
        // duplicate the global min-max pair in 2: one min-saddle pair and
        // one saddle-max pair
        diag[0].death.type = ttk::CriticalType::Saddle1;
        // store the saddle max pair at the vector end
        diag.emplace_back(diag[0]);
        diag.back().birth.type = ttk::CriticalType::Saddle1;
        diag.back().death.type = ttk::CriticalType::Local_maximum;
      }
      max_persistences[i] = diag[0].persistence;
    }

    this->max_dimension_total_
      = *std::max_element(max_persistences.begin(), max_persistences.end());

    if(this->Method == METHOD::PROGRESSIVE) {

      if(!UseInterruptible) {
        TimeLimit = 999999999;
      }

      inv_clustering_ = this->execute(
        intermediateDiagrams_, final_centroids_, all_matchings_);
      needUpdate_ = false;

    } else if(this->Method == METHOD::AUCTION) {

      final_centroids_.resize(1);
      inv_clustering_.resize(numInputs);
      for(int i_input = 0; i_input < numInputs; i_input++) {
        inv_clustering_[i_input] = 0;
      }
      PersistenceDiagramBarycenter pdBarycenter{};

      const auto wassersteinMetric = std::to_string(WassersteinMetric);
      pdBarycenter.setWasserstein(wassersteinMetric);
      pdBarycenter.setMethod(2);
      pdBarycenter.setNumberOfInputs(numInputs);
      pdBarycenter.setTimeLimit(TimeLimit);
      pdBarycenter.setDeterministic(Deterministic);
      pdBarycenter.setUseProgressive(UseProgressive);
      pdBarycenter.setDebugLevel(debugLevel_);
      pdBarycenter.setThreadNumber(threadNumber_);
      pdBarycenter.setAlpha(Alpha);
      pdBarycenter.setLambda(Lambda);

      all_matchings_.resize(1); // at least of size 1
      pdBarycenter.execute(
        intermediateDiagrams_, final_centroids_[0], all_matchings_[0]);

      needUpdate_ = false;
    }
  }

  outputClusteredDiagrams(output_clusters, input, this->intermediateDiagrams_,
                          this->all_matchings_, this->inv_clustering_,
                          this->DisplayMethod, this->Spacing,
                          this->max_dimension_total_);
  outputCentroids(output_centroids, this->final_centroids_,
                  this->all_matchings_, input[0], this->DisplayMethod,
                  this->Spacing, this->max_dimension_total_);
  outputMatchings(
    output_matchings, this->NumberOfClusters, this->intermediateDiagrams_,
    this->all_matchings_, this->final_centroids_, this->inv_clustering_,
    this->DisplayMethod, this->Spacing, this->max_dimension_total_);

  return 1;
}

void addCostsAsFieldData(vtkUnstructuredGrid *vtu,
                         const double minSadCost,
                         const double sadSadCost,
                         const double sadMaxCost) {

  // add global matchings cost as FieldData (only 1 tuple per array)
  vtkNew<vtkDoubleArray> minSad{};
  minSad->SetName("MinSaddleCost");
  minSad->SetNumberOfTuples(1);
  minSad->SetTuple1(0, minSadCost);
  vtu->GetFieldData()->AddArray(minSad);

  vtkNew<vtkDoubleArray> sadSad{};
  sadSad->SetName("SaddleSaddleCost");
  sadSad->SetNumberOfTuples(1);
  sadSad->SetTuple1(0, sadSadCost);
  vtu->GetFieldData()->AddArray(sadSad);

  vtkNew<vtkDoubleArray> sadMax{};
  sadMax->SetName("SaddleMaxCost");
  sadMax->SetNumberOfTuples(1);
  sadMax->SetTuple1(0, sadMaxCost);
  vtu->GetFieldData()->AddArray(sadMax);
}

void ttkPersistenceDiagramClustering::outputClusteredDiagrams(
  vtkMultiBlockDataSet *output,
  const std::vector<vtkUnstructuredGrid *> &diagsVTU,
  const std::vector<ttk::DiagramType> &diags,
  const std::vector<std::vector<std::vector<ttk::MatchingType>>>
    &matchingsPerCluster,
  const std::vector<int> &inv_clustering,
  const DISPLAY dm,
  const double spacing,
  const double max_persistence) const {

  // index of diagram in its cluster
  std::vector<int> diagIdInClust{};
  // number of diagrams per cluster
  std::vector<int> clustSize{};

  // prep work for displaying diagrams as cluster stars
  if(dm == DISPLAY::STARS) {
    // total number of clusters
    const auto nClusters
      = 1 + *std::max_element(inv_clustering.begin(), inv_clustering.end());
    clustSize.resize(nClusters, 0);
    diagIdInClust.resize(diagsVTU.size());
    for(size_t i = 0; i < inv_clustering.size(); ++i) {
      auto &diagsInClust = clustSize[inv_clustering[i]];
      diagIdInClust[i] = diagsInClust;
      diagsInClust++;
    }
  }

  output->SetNumberOfBlocks(diagsVTU.size());

  for(size_t i = 0; i < diagsVTU.size(); ++i) {
    vtkNew<vtkUnstructuredGrid> vtu{};
    vtu->ShallowCopy(diagsVTU[i]);

    vtkNew<vtkIntArray> clusterId{};
    clusterId->SetName("ClusterID");
    clusterId->SetNumberOfComponents(1);
    clusterId->SetNumberOfTuples(vtu->GetNumberOfPoints());
    clusterId->Fill(inv_clustering[i]);
    vtu->GetPointData()->AddArray(clusterId);

    // add clusterId to FieldData too (only 1 tuple)
    vtkNew<vtkIntArray> cidFieldData{};
    cidFieldData->SetName("ClusterID");
    cidFieldData->SetNumberOfComponents(1);
    cidFieldData->SetNumberOfTuples(1);
    cidFieldData->Fill(inv_clustering[i]);
    vtu->GetFieldData()->AddArray(cidFieldData);

    // add Persistence data array on vertices
    vtkNew<vtkDoubleArray> pointPers{};
    pointPers->SetName("Persistence");
    pointPers->SetNumberOfTuples(vtu->GetNumberOfPoints());
    vtu->GetPointData()->AddArray(pointPers);

    // diagonal uses two existing points
    for(int j = 0; j < vtu->GetNumberOfCells() - 1; ++j) {
      const auto persArray = vtu->GetCellData()->GetArray("Persistence");
      const auto pers = persArray->GetTuple1(j);
      pointPers->SetTuple1(2 * j + 0, pers);
      pointPers->SetTuple1(2 * j + 1, pers);
    }

    const auto cid = inv_clustering[i];
    const auto &matchings{matchingsPerCluster[cid][i]};
    double minSadCost{}, sadSadCost{}, sadMaxCost{};
    const auto &diag{diags[i]};

    for(size_t j = 0; j < matchings.size(); ++j) {
      const auto &m{matchings[j]};
      const auto bidderId{std::get<0>(m)};

      // avoid out-of-bound accesses
      if(bidderId >= static_cast<ttk::SimplexId>(diag.size())) {
        this->printWrn("Out-of-bounds access averted");
        continue;
      }

      const auto &p1{diag[bidderId]};
      if(p1.birth.type == ttk::CriticalType::Local_minimum) {
        minSadCost += std::get<2>(m);
      } else if(p1.birth.type == ttk::CriticalType::Saddle1
                && p1.death.type == ttk::CriticalType::Saddle2) {
        sadSadCost += std::get<2>(m);
      } else if(p1.death.type == ttk::CriticalType::Local_maximum) {
        sadMaxCost += std::get<2>(m);
      }
    }

    addCostsAsFieldData(vtu, minSadCost, sadSadCost, sadMaxCost);

    if(dm == DISPLAY::MATCHINGS && spacing > 0) {
      // translate diagrams along the Z axis
      vtkNew<vtkTransform> tr{};
      tr->Translate(0, 0, i == 0 ? -spacing : spacing);
      vtkNew<vtkTransformFilter> trf{};
      trf->SetTransform(tr);
      trf->SetInputData(vtu);
      trf->Update();
      output->SetBlock(i, trf->GetOutputDataObject(0));
    } else if(dm == DISPLAY::STARS && spacing > 0) {
      const auto c = inv_clustering[i];
      const auto angle = 2.0 * M_PI * static_cast<double>(diagIdInClust[i])
                         / static_cast<double>(clustSize[c]);
      // translate diagrams in the XY plane
      vtkNew<vtkTransform> tr{};
      tr->Translate(3.0 * (spacing + 0.2) * max_persistence * c
                      + spacing * max_persistence * std::cos(angle) + 0.2,
                    spacing * max_persistence * std::sin(angle), 0);
      vtkNew<vtkTransformFilter> trf{};
      trf->SetTransform(tr);
      trf->SetInputData(vtu);
      trf->Update();
      output->SetBlock(i, trf->GetOutputDataObject(0));

    } else {
      // add diagram to output multi-block dataset
      output->SetBlock(i, vtu);
    }
  }
}

void ttkPersistenceDiagramClustering::outputCentroids(
  vtkMultiBlockDataSet *output,
  const std::vector<DiagramType> &final_centroids,
  const std::vector<std::vector<std::vector<ttk::MatchingType>>>
    &matchingsPerCluster,
  vtkUnstructuredGrid *const someInputDiag,
  const DISPLAY dm,
  const double spacing,
  const double max_persistence) const {

  if(final_centroids.size() != matchingsPerCluster.size()) {
    this->printWrn("Inconsistent matchings vector size");
  }

  const auto da{
    someInputDiag->GetCellData()->GetArray(ttk::PersistenceBirthName)};
  const auto dim{static_cast<int>(someInputDiag->GetCellData()
                                    ->GetArray(ttk::PersistencePairTypeName)
                                    ->GetRange()[1])
                 + 1};

  for(size_t i = 0; i < final_centroids.size(); ++i) {
    vtkNew<vtkUnstructuredGrid> vtu{};
    DiagramToVTU(vtu, final_centroids[i], da, *this, dim, false);

    vtkNew<vtkIntArray> clusterId{};
    clusterId->SetName("ClusterID");
    clusterId->SetNumberOfTuples(vtu->GetNumberOfPoints());
    clusterId->Fill(i);
    vtu->GetPointData()->AddArray(clusterId);

    // add clusterId to FieldData too (only 1 tuple)
    vtkNew<vtkIntArray> cidFieldData{};
    cidFieldData->SetName("ClusterID");
    cidFieldData->SetNumberOfComponents(1);
    cidFieldData->SetNumberOfTuples(1);
    cidFieldData->Fill(i);
    vtu->GetFieldData()->AddArray(cidFieldData);

    vtkNew<vtkDoubleArray> pointPers{};
    pointPers->SetName("Persistence");
    pointPers->SetNumberOfTuples(vtu->GetNumberOfPoints());
    vtu->GetPointData()->AddArray(pointPers);

    for(size_t j = 0; j < final_centroids[i].size(); ++j) {
      const auto &pair{final_centroids[i][j]};
      pointPers->SetTuple1(2 * j + 0, pair.persistence);
      pointPers->SetTuple1(2 * j + 1, pair.persistence);
    }

    double minSadCost{}, sadSadCost{}, sadMaxCost{};

    for(const auto &matchingsPerDiag : matchingsPerCluster[i]) {
      for(const auto &m : matchingsPerDiag) {
        const auto goodId{std::get<1>(m)};
        const auto &p0{final_centroids[i][goodId]};
        if(p0.birth.type == ttk::CriticalType::Local_minimum) {
          minSadCost += std::get<2>(m);
        } else if(p0.birth.type == ttk::CriticalType::Saddle1
                  && p0.death.type == ttk::CriticalType::Saddle2) {
          sadSadCost += std::get<2>(m);
        } else if(p0.death.type == ttk::CriticalType::Local_maximum) {
          sadMaxCost += std::get<2>(m);
        }
      }
    }

    addCostsAsFieldData(vtu, minSadCost, sadSadCost, sadMaxCost);

    if(dm == DISPLAY::STARS && spacing > 0) {
      // shift centroid along the X axis
      vtkNew<vtkTransform> tr{};
      tr->Translate(3.0 * (spacing + 0.2) * max_persistence * i, 0, 0);
      vtkNew<vtkTransformFilter> trf{};
      trf->SetTransform(tr);
      trf->SetInputData(vtu);
      trf->Update();
      output->SetBlock(i, trf->GetOutputDataObject(0));

    } else {
      // add centroid to output multi-block dataset
      output->SetBlock(i, vtu);
    }
  }
}

void ttkPersistenceDiagramClustering::outputMatchings(
  vtkMultiBlockDataSet *output,
  const size_t nClusters,
  const std::vector<DiagramType> &diags,
  const std::vector<std::vector<std::vector<ttk::MatchingType>>>
    &matchingsPerCluster,
  const std::vector<DiagramType> &centroids,
  const std::vector<int> &inv_clustering,
  const DISPLAY dm,
  const double spacing,
  const double max_persistence) const {

  // index of diagram in its cluster
  std::vector<int> diagIdInClust{};
  // number of diagrams per cluster
  std::vector<int> clustSize{};

  // prep work for displaying diagrams as cluster stars
  if(dm == DISPLAY::STARS) {
    clustSize.resize(nClusters, 0);
    diagIdInClust.resize(diags.size());
    for(size_t i = 0; i < inv_clustering.size(); ++i) {
      auto &diagsInClust = clustSize[inv_clustering[i]];
      diagIdInClust[i] = diagsInClust;
      diagsInClust++;
    }
  }

  // count the number of bidders per centroid pair
  // (when with only 1 cluster and 2 diagrams)
  std::vector<int> matchings_count(centroids[0].size());
  std::vector<int> count_to_good{};

  for(size_t i = 0; i < diags.size(); ++i) {
    const auto cid = inv_clustering[i];
    const auto &diag{diags[i]};
    const auto &matchings{matchingsPerCluster[cid][i]};

    vtkNew<vtkUnstructuredGrid> matchingsGrid{};

    const auto nCells{matchings.size()};
    const auto nPoints{2 * matchings.size()};

    vtkNew<vtkPoints> points{};
    points->SetNumberOfPoints(nPoints);
    matchingsGrid->SetPoints(points);

    // point data
    vtkNew<vtkIntArray> diagIdVerts{};
    diagIdVerts->SetName("DiagramID");
    diagIdVerts->SetNumberOfTuples(nPoints);
    matchingsGrid->GetPointData()->AddArray(diagIdVerts);

    vtkNew<vtkIntArray> pointId{};
    pointId->SetName("PointID");
    pointId->SetNumberOfTuples(nPoints);
    matchingsGrid->GetPointData()->AddArray(pointId);

    // cell data
    vtkNew<vtkIntArray> diagIdCells{};
    diagIdCells->SetName("DiagramID");
    diagIdCells->SetNumberOfTuples(nCells);
    matchingsGrid->GetCellData()->AddArray(diagIdCells);

    vtkNew<vtkIntArray> clusterId{};
    clusterId->SetName("ClusterID");
    clusterId->SetNumberOfTuples(nCells);
    clusterId->Fill(cid);
    matchingsGrid->GetCellData()->AddArray(clusterId);

    vtkNew<vtkDoubleArray> matchCost{};
    matchCost->SetName("Cost");
    matchCost->SetNumberOfTuples(nCells);
    matchingsGrid->GetCellData()->AddArray(matchCost);

    vtkNew<vtkIntArray> pairType{};
    pairType->SetName("PairType");
    pairType->SetNumberOfTuples(nCells);
    matchingsGrid->GetCellData()->AddArray(pairType);

    double minSadCost{}, sadSadCost{}, sadMaxCost{};

    for(size_t j = 0; j < matchings.size(); ++j) {
      const auto &m{matchings[j]};
      const auto bidderId{std::get<0>(m)};
      const auto goodId{std::get<1>(m)};

      // avoid out-of-bound accesses
      if(goodId >= static_cast<ttk::SimplexId>(centroids[cid].size())
         || bidderId >= static_cast<ttk::SimplexId>(diag.size())) {
        this->printWrn("Out-of-bounds access averted");
        continue;
      }

      if(nClusters == 1) {
        matchings_count[goodId] += 1;
        count_to_good.push_back(goodId);
      }

      const auto &p0{centroids[cid][goodId]};
      const auto &p1{diag[bidderId]};
      std::array<double, 3> coords0{p0.birth.sfValue, p0.death.sfValue, 0};
      std::array<double, 3> coords1{p1.birth.sfValue, p1.death.sfValue, 0};

      if(dm == DISPLAY::STARS && spacing > 0) {
        const auto angle = 2.0 * M_PI * static_cast<double>(diagIdInClust[i])
                           / static_cast<double>(clustSize[cid]);
        const auto shift
          = 3.0 * (std::abs(spacing) + 0.2) * max_persistence * cid;
        coords0[0] += shift;
        coords1[0] += shift + spacing * max_persistence * std::cos(angle);
        coords1[1] += spacing * max_persistence * std::sin(angle);

      } else if(dm == DISPLAY::MATCHINGS) {
        coords1[2] = (diags.size() == 2 && i == 0) ? -spacing : spacing;
      }

      points->SetPoint(2 * j + 0, coords0.data());
      points->SetPoint(2 * j + 1, coords1.data());
      std::array<vtkIdType, 2> ids{
        2 * static_cast<vtkIdType>(j) + 0,
        2 * static_cast<vtkIdType>(j) + 1,
      };
      matchingsGrid->InsertNextCell(VTK_LINE, 2, ids.data());

      diagIdCells->SetTuple1(j, i);
      matchCost->SetTuple1(j, std::get<2>(m));
      diagIdVerts->SetTuple1(2 * j + 0, i);
      diagIdVerts->SetTuple1(2 * j + 1, i);
      pointId->SetTuple1(2 * j + 0, goodId);
      pointId->SetTuple1(2 * j + 1, bidderId);
      pairType->SetTuple1(j, p1.dim);

      if(p1.birth.type == ttk::CriticalType::Local_minimum) {
        minSadCost += std::get<2>(m);
      } else if(p1.birth.type == ttk::CriticalType::Saddle1
                && p1.death.type == ttk::CriticalType::Saddle2) {
        sadSadCost += std::get<2>(m);
      } else if(p1.death.type == ttk::CriticalType::Local_maximum) {
        sadMaxCost += std::get<2>(m);
      }
    }

    addCostsAsFieldData(matchingsGrid, minSadCost, sadSadCost, sadMaxCost);

    // add diagram matchings to multi-block
    output->SetBlock(i, matchingsGrid);
  }

  // add matchings number
  if(nClusters == 1 && diags.size() == 2) {
    size_t nPrevCells{};
    for(size_t i = 0; i < diags.size(); ++i) {
      vtkNew<vtkIntArray> matchNumber{};
      matchNumber->SetName("MatchNumber");
      const auto matchings
        = vtkUnstructuredGrid::SafeDownCast(output->GetBlock(i));
      const auto nCells = matchings->GetNumberOfCells();
      matchNumber->SetNumberOfTuples(nCells);
      for(int j = 0; j < nCells; j++) {
        const auto goodId = count_to_good[j + (i == 0 ? 0 : nPrevCells)];
        matchNumber->SetTuple1(j, matchings_count[goodId]);
      }
      matchings->GetCellData()->AddArray(matchNumber);
      nPrevCells = nCells;
    }
  }
}
