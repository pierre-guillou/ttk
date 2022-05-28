#include <Quadrangulation.h>

int ttk::Quadrangulation::preconditionVertexNeighbors() {

  Timer tm;

  this->printMsg(
    "Building vertex neighbors", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  std::vector<std::vector<SimplexId>> vertNeighs(this->nVerts_);

  for(SimplexId i = 0; i < this->nCells_; ++i) {
    const auto &q{this->cells_[i]};
    vertNeighs[q[0]].emplace_back(q[1]);
    vertNeighs[q[0]].emplace_back(q[3]);
    vertNeighs[q[1]].emplace_back(q[0]);
    vertNeighs[q[1]].emplace_back(q[2]);
    vertNeighs[q[2]].emplace_back(q[1]);
    vertNeighs[q[2]].emplace_back(q[3]);
    vertNeighs[q[3]].emplace_back(q[2]);
    vertNeighs[q[3]].emplace_back(q[0]);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < vertNeighs.size(); ++i) {
    auto &vec{vertNeighs[i]};
    std::sort(vec.begin(), vec.end());
    const auto last{std::unique(vec.begin(), vec.end())};
    vec.erase(last, vec.end());
  }

  this->vertexNeighbors_.fillFrom(vertNeighs, this->threadNumber_);

  printMsg("Built " + std::to_string(this->nVerts_) + " vertex neighbors", 1,
           tm.getElapsedTime(), 1);

  return 0;
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

int test() {
  ttk::SurfaceGeometrySmoother surfgeom{};
  ttk::Quadrangulation quad{};
  ttk::Triangulation triangulation{};

  surfgeom.preconditionTriangulationToSmooth(&quad);
  surfgeom.preconditionTriangulationSurface(&triangulation);
  std::vector<float> outputCoords(9);
  std::vector<float> inputCoords(9);

  surfgeom.execute(outputCoords.data(), inputCoords.data(), nullptr, nullptr,
                   10, quad, triangulation);
  return 0;
}
