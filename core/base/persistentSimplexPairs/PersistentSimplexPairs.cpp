#include <PersistentSimplexPairs.h>

ttk::PersistentSimplexPairs::PersistentSimplexPairs() {
  this->setDebugMsgPrefix("PersistentSimplexPairs");
}

template <typename Container>
ttk::SimplexId ttk::PersistentSimplexPairs::eliminateBoundaries(
  const Simplex &c,
  std::vector<bool> &onBoundary,
  std::vector<Container> &boundaries,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder,
  const std::vector<Simplex> &partners) const {

  auto &boundary{boundaries[c.cellId_]};
  this->addCellBoundary(c, onBoundary, boundary);

  const auto getLocalId = [&c, this](const SimplexId a) {
    if(c.dim_ == 1) {
      return a;
    } else if(c.dim_ == 2) {
      return a - this->nVerts_;
    } else if(c.dim_ == 3) {
      return a - this->nVerts_ - this->nEdges_;
    }
    return -1;
  };

  while(!boundary.empty()) {
    // youngest cell on boundary
    const auto tau{getLocalId(*boundary.begin())};
    const Cell cTau{c.dim_ - 1, tau};
    const auto pTau{this->dg_.getPairedCell(cTau)};
    const Simplex *partnerTau{}; // co-facet of tau
    if(pTau == -1) {
      // tau is a critical cell
      partnerTau = &partners[getCellId(c.dim_ - 1, tau)];
    } else {
      partnerTau = &filtration[filtOrder[getCellId(c.dim_, pTau)]];
    }
    if(partnerTau->dim_ == -1 || partnerTau->id_ == -1) {
      return tau;
    }
    if(pTau == -1) {
      for(const auto e : boundaries[partnerTau->cellId_]) {
        if(!onBoundary[e]) {
          boundary.emplace(e);
          onBoundary[e] = true;
        } else {
          const auto it{boundary.find(e)};
          boundary.erase(it);
          onBoundary[e] = false;
        }
      }
    } else {
      this->addCellBoundary(*partnerTau, onBoundary, boundary);
    }
  }

  return -1;
}

int ttk::PersistentSimplexPairs::pairCells(
  std::vector<PersistencePair> &pairs,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder) const {

  // paired simplices
  std::vector<Simplex> partners(filtration.size());
  std::vector<bool> onBoundary(filtration.size(), false);

  const auto cmpSimplices = [&filtOrder](const SimplexId a, const SimplexId b) {
    return filtOrder[a] > filtOrder[b];
  };

  // boundaries storage
  using Container = std::set<SimplexId, decltype(cmpSimplices)>;
  std::vector<Container> boundaries(filtration.size(), Container(cmpSimplices));

  // critical simplices indices in filtration vector
  std::vector<SimplexId> critFilt{};
  for(size_t i = 0; i < filtration.size(); ++i) {
    const auto &s{filtration[i]};
    if(this->dg_.isCellCritical(Cell{s.dim_, s.id_}) && s.dim_ > 0) {
      critFilt.emplace_back(i);
    }
  }

  Timer tm{};

  this->printMsg("Computing pairs", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  for(size_t i = 0; i < critFilt.size(); ++i) {

    const auto &c{filtration[critFilt[i]]};

    const auto partner = eliminateBoundaries(
      c, onBoundary, boundaries, filtration, filtOrder, partners);
    if(partner != -1) {
      const auto &pc{filtration[filtOrder[getCellId(c.dim_ - 1, partner)]]};
      partners[c.cellId_] = pc;
      partners[pc.cellId_] = c;

      // only record pairs with non-null persistence
      if(c.vertsOrder_[0] != pc.vertsOrder_[0]) {
        pairs.emplace_back(partner, c.id_, c.dim_ - 1);
      }
    }

    // clean mask
    for(const auto e : boundaries[c.cellId_]) {
      onBoundary[e] = false;
    }

    if(i % (critFilt.size() / 10) == 0) {
      this->printMsg(
        "Computing pairs",
        std::round(10 * i / static_cast<float>(critFilt.size())) / 10.0f,
        tm.getElapsedTime(), 1, ttk::debug::LineMode::REPLACE);
    }
  }

  const auto nRegPairs{pairs.size()};

  this->printMsg("Computed " + std::to_string(nRegPairs) + " regular pair"
                   + (nRegPairs > 1 ? "s" : ""),
                 1.0, tm.getElapsedTime(), 1);

  // get infinite pairs
  for(SimplexId i = 0; i < this->nVerts_; ++i) {
    if(partners[i].id_ == -1 && this->dg_.isCellCritical(Cell{0, i})) {
      pairs.emplace_back(i, -1, 0);
    }
  }
  if(this->nTri_ > 0) {
    for(SimplexId i = 0; i < this->nEdges_; ++i) {
      if(partners[i + this->nVerts_].id_ == -1
         && this->dg_.isCellCritical(Cell{1, i})) {
        pairs.emplace_back(i, -1, 1);
      }
    }
  }
  if(this->nTetra_ > 0) {
    for(SimplexId i = 0; i < this->nTri_; ++i) {
      if(partners[i + this->nVerts_ + this->nEdges_].id_ == -1
         && this->dg_.isCellCritical(Cell{2, i})) {
        pairs.emplace_back(i, -1, 2);
      }
    }
  }

  const auto nInfPairs{pairs.size() - nRegPairs};

  this->printMsg("Detected " + std::to_string(nInfPairs) + " infinite pair"
                 + (nInfPairs > 1 ? "s" : ""));

  return 0;
}
