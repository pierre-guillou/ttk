#include <QuadrangulationSubdivision.h>

ttk::SimplexId ttk::QuadrangulationSubdivision::findQuadBary(
  const std::vector<size_t> &quadVertices) const {

  std::vector<float> sum(vertexDistance_[*quadVertices.begin()].size(),
                         std::numeric_limits<float>::infinity());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sum.size(); ++i) {

    // skip following computation if too far from any parent quad vertex
    bool skip = false;

    for(const auto vert : quadVertices) {
      if(vertexDistance_[vert][i] == std::numeric_limits<float>::infinity()) {
        skip = true;
        break;
      }
    }

    if(skip) {
      continue;
    }

    float m = vertexDistance_[quadVertices[0]][i];
    float n = vertexDistance_[quadVertices[1]][i];
    float o = vertexDistance_[quadVertices[2]][i];
    float p = vertexDistance_[quadVertices[3]][i];

    // try to be "near" the four parent vertices
    sum[i] = m + n + o + p;

    // try to be on the diagonals intersection
    sum[i] += std::abs(m - o);
    sum[i] += std::abs(n - p);
  }

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

int ttk::QuadrangulationSubdivision::getQuadNeighbors(
  const std::vector<Quad> &quads,
  std::vector<std::set<size_t>> &neighbors,
  const bool secondNeighbors) const {
  Timer tm;

  for(auto &q : quads) {
    auto i = static_cast<size_t>(q[0]);
    auto j = static_cast<size_t>(q[1]);
    auto k = static_cast<size_t>(q[2]);
    auto l = static_cast<size_t>(q[3]);
    if(secondNeighbors) {
      neighbors[i].insert(j);
      neighbors[i].insert(k);
      neighbors[i].insert(l);
      neighbors[j].insert(i);
      neighbors[j].insert(k);
      neighbors[j].insert(l);
      neighbors[k].insert(i);
      neighbors[k].insert(j);
      neighbors[k].insert(l);
      neighbors[l].insert(i);
      neighbors[l].insert(j);
      neighbors[l].insert(k);
    } else {
      neighbors[i].insert(j);
      neighbors[i].insert(l);
      neighbors[k].insert(j);
      neighbors[k].insert(l);
      neighbors[j].insert(k);
      neighbors[j].insert(i);
      neighbors[l].insert(k);
      neighbors[l].insert(i);
    }
  }

  this->printMsg("Computed neighbors mapping of "
                   + std::to_string(outputPoints_.size()) + " points",
                 1.0, tm.getElapsedTime(), debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
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
