#include <PersistentSimplexPairs.h>

ttk::PersistentSimplexPairs::PersistentSimplexPairs() {
  this->setDebugMsgPrefix("PersistentSimplexPairs");
}

ttk::SimplexId ttk::PersistentSimplexPairs::eliminateBoundaries(
  const Simplex &c,
  VisitedMask &boundary,
  const std::vector<SimplexId> &filtOrder,
  const std::vector<Simplex> &partners) const {

  this->addCellBoundary(c, boundary);

  while(!boundary.visitedIds_.empty()) {
    // youngest cell on boundary
    const auto tau{*std::max_element(
      boundary.visitedIds_.begin(), boundary.visitedIds_.end(),
      [&filtOrder, &c, this](const SimplexId a, const SimplexId b) {
        return filtOrder[getCellId(c.dim_ - 1, a)]
               < filtOrder[getCellId(c.dim_ - 1, b)];
      })};
    const auto &partnerTau{partners[getCellId(c.dim_ - 1, tau)]};
    if(partnerTau.dim_ == -1 || partnerTau.id_ == -1) {
      return tau;
    }
    addCellBoundary(partnerTau, boundary);
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
  std::vector<Simplex> partners(filtration.size());

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
      partners[c.cellId_] = pc;
      partners[pc.cellId_] = c;

      // only record pairs with non-null persistence
      if(c.vertsOrder_[0] != pc.vertsOrder_[0]) {
        pairs.emplace_back(partner, c.id_, c.dim_ - 1);
      }
    }

    if(i % (filtration.size() / 10) == 0) {
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
    if(partners[i].id_ == -1) {
      pairs.emplace_back(i, -1, 0);
    }
  }
  if(this->nTri_ > 0) {
    for(SimplexId i = 0; i < this->nEdges_; ++i) {
      if(partners[i + this->nVerts_].id_ == -1) {
        pairs.emplace_back(i, -1, 1);
      }
    }
  }
  if(this->nTetra_ > 0) {
    for(SimplexId i = 0; i < this->nTri_; ++i) {
      if(partners[i + this->nVerts_ + this->nEdges_].id_ == -1) {
        pairs.emplace_back(i, -1, 2);
      }
    }
  }

  const auto nInfPairs{pairs.size() - nRegPairs};

  this->printMsg("Detected " + std::to_string(nInfPairs) + " infinite pair"
                 + (nInfPairs > 1 ? "s" : ""));

  return 0;
}

template <typename Container>
ttk::SimplexId ttk::PersistentSimplexPairs::eliminateBoundariesV4(
  const Simplex &c,
  std::vector<bool> &onBoundary,
  std::vector<Container> &boundaries,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder,
  const std::vector<Simplex> &partners) const {

  auto &boundary{boundaries[c.cellId_]};
  if(!boundary.empty()) {
    for(const auto e : boundary) {
      onBoundary[e] = true;
    }
  } else {
    this->addCellBoundary(c, onBoundary, boundary);
  }

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

int ttk::PersistentSimplexPairs::pairCellsV4(
  std::vector<PersistencePair> &pairs,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder) const {

  Timer tmall{};

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
  std::array<std::vector<SimplexId>, 3> critFilt{};
  for(size_t i = 0; i < filtration.size(); ++i) {
    const auto &s{filtration[i]};
    if(this->dg_.isCellCritical(Cell{s.dim_, s.id_}) && s.dim_ > 0) {
      critFilt[s.dim_ - 1].emplace_back(i);
    }
  }

  this->printMsg("Memory allocations", 1.0, tmall.getElapsedTime(), 1,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  Timer tmpar{};

  for(const auto &vec : critFilt) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic) \
  firstprivate(onBoundary)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < vec.size(); ++i) {
      const auto &c{filtration[vec[i]]};

      this->eliminateBoundariesV4(
        c, onBoundary, boundaries, filtration, filtOrder, partners);

      // clean mask
      for(const auto e : boundaries[c.cellId_]) {
        onBoundary[e] = false;
      }
    }
  }

  Timer tm{};

  const auto processDim = [&](const std::vector<SimplexId> &critSimplices,
                              std::vector<PersistencePair> &res) {
    for(size_t i = 0; i < critSimplices.size(); ++i) {

      const auto &c{filtration[critSimplices[i]]};

      const auto partner = eliminateBoundariesV4(
        c, onBoundary, boundaries, filtration, filtOrder, partners);
      if(partner != -1) {
        const auto &pc{filtration[filtOrder[getCellId(c.dim_ - 1, partner)]]};
        partners[c.cellId_] = pc;
        partners[pc.cellId_] = c;

        // only record pairs with non-null persistence
        if(c.vertsOrder_[0] != pc.vertsOrder_[0]) {
          res.emplace_back(partner, c.id_, c.dim_ - 1);
        }
      }

      // clean mask
      for(const auto e : boundaries[c.cellId_]) {
        onBoundary[e] = false;
      }
    }
  };

  const auto dim{this->dg_.getDimensionality()};

  // avoid concurrent writes into pairs vector
  std::vector<PersistencePair> pairsSaddleMax{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections if(dim > 1)
#endif // TTK_ENABLE_OPENMP
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif // TTK_ENABLE_OPENMP
    {
      Timer tcrit{};
      processDim(critFilt[0], pairs);
      this->printMsg(
        "Computed min-saddle pairs", 1.0, tcrit.getElapsedTime(), 1);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif // TTK_ENABLE_OPENMP
    if(dim > 1) {
      Timer tcrit{};
      processDim(critFilt[dim - 1], pairsSaddleMax);
      this->printMsg(
        "Computed saddle-max pairs", 1.0, tcrit.getElapsedTime(), 1);
    }
  }

  pairs.insert(pairs.end(), pairsSaddleMax.begin(), pairsSaddleMax.end());

  if(dim > 2) { // sandwich
    Timer tcrit{};
    std::vector<SimplexId> nonPaired2Saddles{};
    nonPaired2Saddles.reserve(critFilt[1].size());
    for(const auto s2 : critFilt[1]) {
      const auto &p{partners[filtration[s2].cellId_]};
      if(p.id_ == -1 || p.dim_ == -1) {
        nonPaired2Saddles.emplace_back(s2);
      }
    }
    processDim(nonPaired2Saddles, pairs);
    this->printMsg(
      "Computed saddle-saddle pairs", 1.0, tcrit.getElapsedTime(), 1);
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

template <typename Container>
ttk::SimplexId ttk::PersistentSimplexPairs::eliminateBoundariesV3(
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

int ttk::PersistentSimplexPairs::pairCellsV3(
  std::vector<PersistencePair> &pairs,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder) const {

  Timer tmall{};

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
  std::array<std::vector<SimplexId>, 3> critFilt{};
  for(size_t i = 0; i < filtration.size(); ++i) {
    const auto &s{filtration[i]};
    if(this->dg_.isCellCritical(Cell{s.dim_, s.id_}) && s.dim_ > 0) {
      critFilt[s.dim_ - 1].emplace_back(i);
    }
  }

  this->printMsg("Memory allocations", 1.0, tmall.getElapsedTime(), 1,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  Timer tm{};

  const auto processDim = [&](const std::vector<SimplexId> &critSimplices) {
    for(size_t i = 0; i < critSimplices.size(); ++i) {

      const auto &c{filtration[critSimplices[i]]};

      const auto partner = eliminateBoundariesV3(
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
    }
  };

  const auto dim{this->dg_.getDimensionality()};

  {
    Timer tcrit{};
    processDim(critFilt[0]);
    this->printMsg("Computed min-saddle pairs", 1.0, tcrit.getElapsedTime(), 1);
  }

  if(dim > 1) {
    Timer tcrit{};
    processDim(critFilt[dim - 1]);
    this->printMsg("Computed saddle-max pairs", 1.0, tcrit.getElapsedTime(), 1);
  }

  if(dim > 2) { // sandwich
    Timer tcrit{};
    std::vector<SimplexId> nonPaired2Saddles{};
    nonPaired2Saddles.reserve(critFilt[1].size());
    for(const auto s2 : critFilt[1]) {
      const auto &p{partners[filtration[s2].cellId_]};
      if(p.id_ == -1 || p.dim_ == -1) {
        nonPaired2Saddles.emplace_back(s2);
      }
    }
    processDim(nonPaired2Saddles);
    this->printMsg(
      "Computed saddle-saddle pairs", 1.0, tcrit.getElapsedTime(), 1);
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

template <typename Container>
ttk::SimplexId ttk::PersistentSimplexPairs::eliminateBoundariesV2(
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

int ttk::PersistentSimplexPairs::pairCellsV2(
  std::vector<PersistencePair> &pairs,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder) const {

  Timer tmall{};

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

  this->printMsg("Memory allocations", 1.0, tmall.getElapsedTime(), 1,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  Timer tm{};

  this->printMsg("Computing pairs", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  for(size_t i = 0; i < critFilt.size(); ++i) {

    const auto &c{filtration[critFilt[i]]};

    const auto partner = eliminateBoundariesV2(
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

ttk::SimplexId ttk::PersistentSimplexPairs::eliminateBoundariesV1(
  const Simplex &c,
  VisitedMask &boundary,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder,
  const std::vector<Simplex> &partners) const {

  this->addCellBoundary(c, boundary);

  while(!boundary.visitedIds_.empty()) {
    // youngest cell on boundary
    const auto tau{*std::max_element(
      boundary.visitedIds_.begin(), boundary.visitedIds_.end(),
      [&filtOrder, &c, this](const SimplexId a, const SimplexId b) {
        return filtOrder[getCellId(c.dim_ - 1, a)]
               < filtOrder[getCellId(c.dim_ - 1, b)];
      })};
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
    addCellBoundary(*partnerTau, boundary);
  }

  return -1;
}

int ttk::PersistentSimplexPairs::pairCellsV1(
  std::vector<PersistencePair> &pairs,
  std::array<std::vector<bool>, 3> &boundaries,
  const std::vector<Simplex> &filtration,
  const std::vector<SimplexId> &filtOrder) const {

  // for VisitedMask
  std::vector<SimplexId> visitedIds{};
  // paired simplices
  std::vector<Simplex> partners(filtration.size());

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

    // store the boundary cells
    VisitedMask vm{boundaries[c.dim_ - 1], visitedIds};

    const auto partner
      = eliminateBoundariesV1(c, vm, filtration, filtOrder, partners);
    if(partner != -1) {
      const auto &pc{filtration[filtOrder[getCellId(c.dim_ - 1, partner)]]};
      partners[c.cellId_] = pc;
      partners[pc.cellId_] = c;

      // only record pairs with non-null persistence
      if(c.vertsOrder_[0] != pc.vertsOrder_[0]) {
        pairs.emplace_back(partner, c.id_, c.dim_ - 1);
      }
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
