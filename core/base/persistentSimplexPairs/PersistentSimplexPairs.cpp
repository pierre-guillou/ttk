#include <PersistentSimplexPairs.h>

ttk::PersistentSimplexPairs::PersistentSimplexPairs() {
  this->setDebugMsgPrefix("PersistentSimplexPairs");
}

ttk::SimplexId ttk::PersistentSimplexPairs::eliminateBoundaries(
  const Simplex &c,
  VisitedMask &boundary,
  const std::vector<SimplexId> &filtOrder,
  const std::vector<const Simplex *> &partners) const {

  this->addCellBoundary(c, boundary);

  while(!boundary.visitedIds_.empty()) {
    // youngest cell on boundary
    const auto tau{*std::max_element(
      boundary.visitedIds_.begin(), boundary.visitedIds_.end(),
      [&filtOrder, &c, this](const SimplexId a, const SimplexId b) {
        return filtOrder[getCellId(c.dim_ - 1, a)]
               < filtOrder[getCellId(c.dim_ - 1, b)];
      })};
    const auto partnerTau{partners[getCellId(c.dim_ - 1, tau)]};
    if(partnerTau == nullptr) {
      return tau;
    }
    addCellBoundary(*partnerTau, boundary);
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
  std::vector<const Simplex *> partners(filtration.size());

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

    const auto partner = eliminateBoundaries(c, vm, filtOrder, partners);
    if(partner != -1) {
      const auto &pc{filtration[filtOrder[getCellId(c.dim_ - 1, partner)]]};
      partners[c.cellId_] = &pc;
      partners[pc.cellId_] = &c;
    }

    if(filtration.size() > 10 && i % (filtration.size() / 10) == 0) {
      this->printMsg(
        "Computing pairs",
        std::round(10 * i / static_cast<float>(filtration.size())) / 10.0f,
        tm.getElapsedTime(), 1, ttk::debug::LineMode::REPLACE);
    }
  }

  std::vector<bool> paired(partners.size(), false);

  for(size_t i = 0; i < partners.size(); ++i) {
    if(paired[i] || partners[i] == nullptr) {
      continue;
    }
    const auto &pc{*partners[i]};
    const auto &c{*partners[pc.cellId_]};

    // skill zero-persistence pairs
    if(c.vertsOrder_[0] == pc.vertsOrder_[0]) {
      continue;
    }

    if(pc.dim_ < c.dim_) {
      pairs.emplace_back(pc.id_, c.id_, pc.dim_);
    } else {
      pairs.emplace_back(c.id_, pc.id_, c.dim_);
    }

    paired[i] = true;
    paired[pc.cellId_] = true;
  }

  const auto nRegPairs{pairs.size()};

  this->printMsg("Computed " + std::to_string(nRegPairs) + " regular pair"
                   + (nRegPairs > 1 ? "s" : ""),
                 1.0, tm.getElapsedTime(), 1);

  // get infinite pairs
  for(SimplexId i = 0; i < this->nVerts_; ++i) {
    if(partners[i] == nullptr) {
      pairs.emplace_back(i, -1, 0);
    }
  }
  if(this->nTri_ > 0) {
    for(SimplexId i = 0; i < this->nEdges_; ++i) {
      if(partners[i + this->nVerts_] == nullptr) {
        pairs.emplace_back(i, -1, 1);
      }
    }
  }
  if(this->nTetra_ > 0) {
    for(SimplexId i = 0; i < this->nTri_; ++i) {
      if(partners[i + this->nVerts_ + this->nEdges_] == nullptr) {
        pairs.emplace_back(i, -1, 2);
      }
    }
  }

  const auto nInfPairs{pairs.size() - nRegPairs};

  this->printMsg("Detected " + std::to_string(nInfPairs) + " infinite pair"
                 + (nInfPairs > 1 ? "s" : ""));

  return 0;
}
