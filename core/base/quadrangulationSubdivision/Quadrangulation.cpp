#include <Quadrangulation.h>

#include <OneSkeleton.h>

ttk::Quadrangulation::Quadrangulation() {
  this->setDebugMsgPrefix("Quadrangulation");
}

int ttk::Quadrangulation::preconditionVertexNeighbors() {

  Timer tm;

  this->preconditionEdges();

  this->printMsg(
    "Building vertex neighbors", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  std::vector<SimplexId> offsets(this->nVerts_ + 1);
  // number of neighbors processed per vertex
  std::vector<SimplexId> neighborsId(this->nVerts_);

  // store number of neighbors per vertex
  for(const auto &e : this->edges_) {
    offsets[e[0] + 1]++;
    offsets[e[1] + 1]++;
  }

  // compute partial sum of number of neighbors per vertex
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat neighbors vector
  std::vector<SimplexId> neighbors(offsets.back());

  // fill flat neighbors vector using offsets and neighbors count vectors
  for(const auto &e : this->edges_) {
    neighbors[offsets[e[0]] + neighborsId[e[0]]] = e[1];
    neighborsId[e[0]]++;
    neighbors[offsets[e[1]] + neighborsId[e[1]]] = e[0];
    neighborsId[e[1]]++;
  }

  // fill FlatJaggedArray struct
  this->vertexNeighbors_.setData(std::move(neighbors), std::move(offsets));

  printMsg("Built " + std::to_string(this->nVerts_) + " vertex neighbors", 1,
           tm.getElapsedTime(), 1);

  return 0;
}

int ttk::Quadrangulation::preconditionVertexStars() {
  Timer tm{};

  printMsg("Building vertex stars", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  std::vector<SimplexId> offsets(this->nVerts_ + 1);
  // number of cells processed per vertex
  std::vector<SimplexId> cellIds(this->nVerts_);

  const auto cellNumber{this->nCells_};

  // store number of stars per vertex
  for(SimplexId i = 0; i < this->nCells_; ++i) {
    const auto &q{this->cells_[i]};
    for(const auto &v : q) {
      offsets[v + 1]++;
    }
  }

  // compute partial sum of number of stars per vertex
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat data vector
  std::vector<SimplexId> data(offsets.back());

  // fill flat data vector using offsets and edges count vectors
  for(SimplexId i = 0; i < cellNumber; ++i) {
    const auto &q{this->cells_[i]};
    for(const auto v : q) {
      data[offsets[v] + cellIds[v]] = i;
      cellIds[v]++;
    }
  }

  // fill FlatJaggedArray struct
  this->vertexStars_.setData(std::move(data), std::move(offsets));

  this->printMsg("Built " + std::to_string(this->nVerts_) + " vertex stars", 1,
                 tm.getElapsedTime(), 1);

  return 0;
}

int ttk::Quadrangulation::preconditionEdges() {

  std::vector<LongSimplexId> offsets(this->nCells_ + 1);
  offsets[0] = 0;
  for(SimplexId i = 0; i < this->nCells_; ++i) {
    offsets[i + 1] = this->cells_[i].size() * (i + 1);
  }
  CellArray ca{this->cells_[0].data(), offsets.data(),
               static_cast<LongSimplexId>(this->nCells_)};
  OneSkeleton osk{};
  osk.setDebugLevel(this->debugLevel_);
  osk.setThreadNumber(this->threadNumber_);

  return osk.buildEdgeList(
    this->nVerts_, ca, this->edges_, this->edgeStars_, this->quadEdges_);
}

void ttk::Quadrangulation::computeStatistics(
  std::vector<SimplexId> &vertsValence,
  std::vector<float> &quadArea,
  std::vector<float> &quadDiagsRatio,
  std::vector<float> &quadEdgesRatio,
  std::vector<float> &quadAnglesRatio) const {

  Timer tm;

  vertsValence.resize(this->nVerts_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nVerts_; ++i) {
    vertsValence[i] = this->getVertexNeighborNumber(i);
  }

  quadArea.resize(this->nCells_);
  quadDiagsRatio.resize(this->nCells_);
  quadEdgesRatio.resize(this->nCells_);
  quadAnglesRatio.resize(this->nCells_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nCells_; ++i) {
    const auto &q = this->cells_[i];
    const auto &pi = this->vertCoords_[q[0]];
    const auto &pj = this->vertCoords_[q[1]];
    const auto &pk = this->vertCoords_[q[2]];
    const auto &pl = this->vertCoords_[q[3]];

    // quadrangle area
    float area0{}, area1{};
    Geometry::computeTriangleArea(pi.data(), pj.data(), pk.data(), area0);
    Geometry::computeTriangleArea(pi.data(), pk.data(), pl.data(), area1);
    quadArea[i] = area0 + area1;

    // diagonals ratio
    const auto diag0 = Geometry::distance(pi.data(), pk.data());
    const auto diag1 = Geometry::distance(pj.data(), pl.data());
    quadDiagsRatio[i] = std::min(diag0, diag1) / std::max(diag0, diag1);

    // edges ratio
    const std::array<float, 4> edges{
      Geometry::distance(pi.data(), pj.data()), // ij
      Geometry::distance(pj.data(), pk.data()), // jk
      Geometry::distance(pk.data(), pl.data()), // kl
      Geometry::distance(pl.data(), pi.data()), // li
    };
    quadEdgesRatio[i] = *std::min_element(edges.begin(), edges.end())
                        / *std::max_element(edges.begin(), edges.end());

    // angles ratio
    const std::array<float, 4> angles{
      Geometry::angle(pi.data(), pl.data(), pi.data(), pj.data()), // lij
      Geometry::angle(pj.data(), pi.data(), pj.data(), pk.data()), // ijk
      Geometry::angle(pk.data(), pj.data(), pk.data(), pl.data()), // jkl
      Geometry::angle(pl.data(), pk.data(), pl.data(), pi.data()), // kli
    };

    const auto min_max{std::minmax_element(angles.begin(), angles.end())};
    quadAnglesRatio[i] = *min_max.first / *min_max.second;
  }

  // compute ratio between quad area and mean quad area

  // global surface area
  float sumArea{};
  for(const auto a : quadArea) {
    sumArea += a;
  }
  for(auto &a : quadArea) {
    a *= quadArea.size() / sumArea;
  }

  this->printMsg("Computed statistics", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}
