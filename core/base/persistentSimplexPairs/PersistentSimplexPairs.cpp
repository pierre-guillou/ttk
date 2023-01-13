#include <PersistentSimplexPairs.h>

ttk::PersistentSimplexPairs::PersistentSimplexPairs() {
  this->setDebugMsgPrefix("PersistentSimplexPairs");
}

ttk::SimplexId ttk::PersistentSimplexPairs::eliminateBoundaries(
  const Simplex &c,
  VisitedMask &boundary,
  std::vector<SimplexId> &partners,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder) const {

  this->addCellBoundary(c, boundary);

  while(!boundary.visitedIds_.empty()) {
    // youngest cell on boundary
    const auto tau{*std::max_element(
      boundary.visitedIds_.begin(), boundary.visitedIds_.end(),
      [&filtOrder, &c, this](const SimplexId a, const SimplexId b) {
        return filtOrder[getCellId(c.dim_ - 1, a)]
               < filtOrder[getCellId(c.dim_ - 1, b)];
      })};
    const auto cTau{getCellId(c.dim_ - 1, tau)};
    const auto partnerTau{partners[cTau]};
    if(partnerTau == -1) {
      partners[c.cellId_] = cTau;
      partners[cTau] = c.cellId_;
      return tau;
    }
    this->addCellBoundary(filtration[filtOrder[partnerTau]], boundary);
  }

  return -1;
}

int ttk::PersistentSimplexPairs::pairCells(
  std::vector<PersistencePair> &pairs,
  std::array<std::vector<bool>, 3> &boundaries,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder) const {

  // for VisitedMask
  std::vector<SimplexId> visitedIds{};
  // paired simplices
  std::vector<SimplexId> partners(filtration.size(), -1);

  Timer tm{};

  this->printMsg("Computing pairs", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  for(size_t i = 0; i < filtration.size(); ++i) {

    const auto &c{filtration[i]};

    // skip vertices
    if(c.dim_ == 0) {
      continue;
    }

    // store the boundary cells
    VisitedMask vm{boundaries[c.dim_ - 1], visitedIds};
    const auto tau
      = this->eliminateBoundaries(c, vm, partners, filtration, filtOrder);
    if(tau != -1) {
      const auto cTau{getCellId(c.dim_ - 1, tau)};
      const auto pc{filtration[filtOrder[cTau]]};
      // only record pairs with non-null persistence
      if(c.vertsOrder_[0] != pc.vertsOrder_[0]) {
        pairs.emplace_back(tau, c.id_, c.dim_ - 1);
      }
    }

    if(filtration.size() > 10 && i % (filtration.size() / 10) == 0) {
      this->printMsg(
        "Computing pairs",
        std::round(10 * i / static_cast<float>(filtration.size())) / 10.0f,
        tm.getElapsedTime(), 1, ttk::debug::LineMode::REPLACE);
    }
  }

  const auto nRegPairs{pairs.size()};

  this->printMsg("Computed " + std::to_string(nRegPairs) + " regular pair"
                   + (nRegPairs > 1 ? "s" : ""),
                 1.0, tm.getElapsedTime(), 1);

  // get infinite pairs
  for(SimplexId i = 0; i < this->nVerts_; ++i) {
    if(partners[i] == -1) {
      pairs.emplace_back(i, -1, 0);
    }
  }
  if(this->nTri_ > 0) {
    for(SimplexId i = 0; i < this->nEdges_; ++i) {
      if(partners[i + this->nVerts_] == -1) {
        pairs.emplace_back(i, -1, 1);
      }
    }
  }
  if(this->nTetra_ > 0) {
    for(SimplexId i = 0; i < this->nTri_; ++i) {
      if(partners[i + this->nVerts_ + this->nEdges_] == -1) {
        pairs.emplace_back(i, -1, 2);
      }
    }
  }

  const auto nInfPairs{pairs.size() - nRegPairs};

  this->printMsg("Detected " + std::to_string(nInfPairs) + " infinite pair"
                 + (nInfPairs > 1 ? "s" : ""));

  return 0;
}
