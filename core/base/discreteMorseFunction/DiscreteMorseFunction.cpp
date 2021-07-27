#include <DiscreteMorseFunction.h>

std::vector<ttk::SimplexId> ttk::DiscreteMorseFunction::topologicalSort(
  const ttk::DiscreteMorseFunction::MonoGraph &monoGraph) const {

  // Kahn algorithm implementation
  std::vector<SimplexId> sortedVertices{};

  // seeds: local maxima (no edges pointing to them)
  std::queue<SimplexId> seeds{};
  const auto nVerts = monoGraph.size();

  std::vector<size_t> nPointed(nVerts);
  for(const auto &vec : monoGraph) {
    for(const auto v : vec) {
      nPointed[v]++;
    }
  }
  for(size_t i = 0; i < nVerts; ++i) {
    if(nPointed[i] == 0) {
      seeds.push(i);
    }
  }

  while(!seeds.empty()) {
    const auto curr = seeds.front();
    seeds.pop();
    sortedVertices.emplace_back(curr);
    for(const auto m : monoGraph[curr]) {
      nPointed[m]--;
      if(nPointed[m] == 0) {
        seeds.push(m);
      }
    }
  }

  if(std::any_of(nPointed.begin(), nPointed.end(),
                 [](const size_t a) { return a > 0; })) {
    this->printErr("Topological sort: cycle in graph");
    return {};
  }

  return sortedVertices;
}
