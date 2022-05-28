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
