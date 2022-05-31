#include <QuadrangulationSubdivision.h>

#include <OneSkeleton.h>

ttk::SimplexId
  ttk::QuadrangulationSubdivision::findQuadBary(const Quad &quad) const {

  std::vector<float> sum(
    vertexDistance_[quad[0]].size(), std::numeric_limits<float>::infinity());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sum.size(); ++i) {

    // skip following computation if too far from any parent quad vertex
    bool skip = false;

    for(const auto vert : quad) {
      if(vertexDistance_[vert][i] == std::numeric_limits<float>::infinity()) {
        skip = true;
        break;
      }
    }

    if(skip) {
      continue;
    }

    float m = vertexDistance_[quad[0]][i];
    float n = vertexDistance_[quad[1]][i];
    float o = vertexDistance_[quad[2]][i];
    float p = vertexDistance_[quad[3]][i];

    // try to be "near" the four parent vertices
    sum[i] = m + n + o + p;

    // try to be on the diagonals intersection
    sum[i] += std::abs(m - o);
    sum[i] += std::abs(n - p);
  }

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

int ttk::QuadrangulationSubdivision::getQuadExtNeighbors(
  FlatJaggedArray &extNeighbors, const std::vector<Quad> &quads) const {

  Timer tm{};
  std::vector<std::vector<SimplexId>> neighbors(this->vertexNumber_);

  for(auto &q : quads) {
    auto i = q[0];
    auto j = q[1];
    auto k = q[2];
    auto l = q[3];
    neighbors[i].emplace_back(j);
    neighbors[i].emplace_back(k);
    neighbors[i].emplace_back(l);
    neighbors[j].emplace_back(i);
    neighbors[j].emplace_back(k);
    neighbors[j].emplace_back(l);
    neighbors[k].emplace_back(i);
    neighbors[k].emplace_back(j);
    neighbors[k].emplace_back(l);
    neighbors[l].emplace_back(i);
    neighbors[l].emplace_back(j);
    neighbors[l].emplace_back(k);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < neighbors.size(); ++i) {
    auto &vec{neighbors[i]};
    std::sort(vec.begin(), vec.end());
    const auto last{std::unique(vec.begin(), vec.end())};
    vec.erase(last, vec.end());
  }

  extNeighbors.fillFrom(neighbors, this->threadNumber_);

  this->printMsg("Computed neighbors mapping of "
                   + std::to_string(outputPoints_.size()) + " points",
                 1.0, tm.getElapsedTime(), debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
}

int ttk::QuadrangulationSubdivision::buildEdgeList(
  std::vector<std::array<SimplexId, 2>> &edgeList,
  FlatJaggedArray &edgeStars,
  std::vector<std::array<SimplexId, 4>> &quadEdgeList,
  const std::vector<Quad> &quads,
  const std::vector<Point> &points) const {

  std::vector<LongSimplexId> offsets(quads.size() + 1);
  offsets[0] = 0;
  for(size_t i = 0; i < quads.size(); ++i) {
    offsets[i + 1] = quads[i].size() * (i + 1);
  }
  CellArray ca{
    quads[0].data(), offsets.data(), static_cast<LongSimplexId>(quads.size())};
  OneSkeleton osk{};
  osk.setDebugLevel(this->debugLevel_);
  osk.setThreadNumber(this->threadNumber_);

  return osk.buildEdgeList(
    points.size(), ca, edgeList, edgeStars, quadEdgeList);
}

void ttk::QuadrangulationSubdivision::clearData() {
  outputQuads_.clear();
  outputPoints_.clear();
  outputValences_.clear();
  outputVertType_.clear();
  outputSubdivision_.clear();
  quadNeighbors_.clear();
  vertexDistance_.clear();
}
