/// \ingroup baseCode
/// \class ttk::dcg::DiscreteMorseFunction
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Yizhe Wang <wangyizhe3518@gmail.com>
/// \date November 2016.
///
/// \brief TTK %discreteMorseFunction processing package.
///
/// %DiscreteMorseFunction is a TTK processing package that handles discrete
/// gradient (in the sense of Discrete Morse Theory).
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <DiscreteGradient.h>

namespace ttk {

  class DiscreteMorseFunction : virtual public dcg::DiscreteGradient {

  public:
    DiscreteMorseFunction() {
      this->setDebugMsgPrefix("DiscreteMorseFunction");
    }

    /**
     * @brief Compute the Discrete Morse Function from the Discrete
     * Gradient
     *
     * @param[out] dmf Discrete Morse Function output: one order
     * value per simplex in the simplicial complex (lowest values:
     * minima, highest values: maxima)
     * @param[in] triangulation Explicit triangulation of the
     * simplicial complex
     *
     * @return 0 in case of success, negative value in case of error
     * @pre DiscreteMorseFunction::buildGradient() needs to be called
     * prior to this function
     */
    template <typename triangulationType>
    std::vector<SimplexId> computeDiscreteMorseFunction(
      const SimplexId *const offset,
      const triangulationType &triangulation) const;

    struct Simplex;

    template <typename triangulationType>
    std::vector<Simplex>
      computeFiltrationOrder(const SimplexId *const offset,
                             const triangulationType &triangulation) const;

    struct Simplex {
      SimplexId dim_{}; // dimension
      SimplexId id_{}; // id in triangulation (overlap between dimensions)
      SimplexId cellId_{}; // cell id (unique)
      // order on vertices, sorted in descending order
      std::array<SimplexId, 4> vertsOrder_{-1, -1, -1, -1};

      friend std::ostream &operator<<(std::ostream &os, const Simplex &rhs) {
        os << rhs.dim_ << " " << rhs.id_ << " (";
        for(SimplexId i = 0; i < rhs.dim_; ++i) {
          os << rhs.vertsOrder_[i] << ", ";
        }
        os << rhs.vertsOrder_[rhs.dim_] << ")";
        return os;
      }

      friend bool operator<(const Simplex &lhs, const Simplex &rhs) {
        return lhs.vertsOrder_ < rhs.vertsOrder_;
      }

      void fillVert(const SimplexId v, const SimplexId *const offset) {
        this->dim_ = 0;
        this->id_ = v;
        this->cellId_ = v;
        this->vertsOrder_[0] = offset[v];
      }

      template <typename triangulationType>
      void fillEdge(const SimplexId e,
                    const SimplexId c,
                    const SimplexId *const offset,
                    const triangulationType &triangulation) {
        this->dim_ = 1;
        this->id_ = e;
        this->cellId_ = c;
        triangulation.getEdgeVertex(e, 0, this->vertsOrder_[0]);
        triangulation.getEdgeVertex(e, 1, this->vertsOrder_[1]);
        this->vertsOrder_[0] = offset[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offset[this->vertsOrder_[1]];
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }

      template <typename triangulationType>
      void fillTriangle(const SimplexId t,
                        const SimplexId c,
                        const SimplexId *const offset,
                        const triangulationType &triangulation) {
        this->dim_ = 2;
        this->id_ = t;
        this->cellId_ = c;
        triangulation.getTriangleVertex(t, 0, this->vertsOrder_[0]);
        triangulation.getTriangleVertex(t, 1, this->vertsOrder_[1]);
        triangulation.getTriangleVertex(t, 2, this->vertsOrder_[2]);
        this->vertsOrder_[0] = offset[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offset[this->vertsOrder_[1]];
        this->vertsOrder_[2] = offset[this->vertsOrder_[2]];
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }

      template <typename triangulationType>
      void fillTetra(const SimplexId T,
                     const SimplexId c,
                     const SimplexId *const offset,
                     const triangulationType &triangulation) {
        this->dim_ = 3;
        this->id_ = T;
        this->cellId_ = c;
        triangulation.getCellVertex(T, 0, this->vertsOrder_[0]);
        triangulation.getCellVertex(T, 1, this->vertsOrder_[1]);
        triangulation.getCellVertex(T, 2, this->vertsOrder_[2]);
        triangulation.getCellVertex(T, 3, this->vertsOrder_[3]);
        this->vertsOrder_[0] = offset[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offset[this->vertsOrder_[1]];
        this->vertsOrder_[2] = offset[this->vertsOrder_[2]];
        this->vertsOrder_[3] = offset[this->vertsOrder_[3]];
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }
    };

  private:
    /**
     * Monotonicity graph: for each vertex, store the edges pointing to
     * the lower neighbors.
     */
    using MonoGraph = std::vector<std::vector<SimplexId>>;

    std::vector<ttk::SimplexId>
      topologicalSort(const MonoGraph &monoGraph) const;
  };
} // namespace ttk

template <typename triangulationType>
std::vector<SimplexId> ttk::DiscreteMorseFunction::computeDiscreteMorseFunction(
  const SimplexId *const offset, const triangulationType &triangulation) const {

  Timer tm{};
  std::vector<SimplexId> res{};

  if(this->gradient_[0].empty()) {
    this->printErr("Empty gradient");
    return res;
  }

  const auto dim = triangulation.getDimensionality();
  const auto nVerts = triangulation.getNumberOfVertices();
  const auto nEdges = triangulation.getNumberOfEdges();
  const auto nTri = dim > 1 ? triangulation.getNumberOfTriangles() : 0;
  const auto nTetra = dim > 2 ? triangulation.getNumberOfCells() : 0;

  // global number of cells (vertices + edges + triangles + tetras)
  const auto num_simplices = nVerts + nEdges + nTri + nTetra;
  // monotonicity graph: each simplex points to its paired co-facet in
  // the Discrete Gradient
  MonoGraph mg(num_simplices);

  // vertices sorted by the inputOffsets field
  std::vector<SimplexId> orderVerts(nVerts);
  for(size_t i = 0; i < orderVerts.size(); ++i) {
    orderVerts[offset[i]] = i;
  }

  // compare consecutive sorted vertices
  for(int i = 0; i < nVerts - 1; ++i) {
    mg[orderVerts[i + 1]].emplace_back(orderVerts[i]);
  }

  for(int i = 0; i < nVerts; ++i) {
    const auto pe = this->getPairedCell(Cell{0, i}, triangulation);
    if(pe != -1) {
      mg[i].emplace_back(nVerts + pe);
    }
    // every other neighbor edge is higher that the current vertex
    const auto ne = triangulation.getVertexEdgeNumber(i);
    for(SimplexId j = 0; j < ne; ++j) {
      SimplexId e{};
      triangulation.getVertexEdge(i, j, e);
      // skip paired edge
      if(e == pe) {
        continue;
      }
      mg[nVerts + e].emplace_back(i);
    }
  }
  for(int i = 0; i < nEdges; ++i) {
    const auto pt = this->getPairedCell(Cell{1, i}, triangulation);
    if(pt != -1) {
      mg[nVerts + i].emplace_back(nVerts + nEdges + pt);
    }
    // every other neighbor triangle is higher that the current edge
    const auto nt = triangulation.getEdgeTriangleNumber(i);
    for(SimplexId j = 0; j < nt; ++j) {
      SimplexId t{};
      triangulation.getEdgeTriangle(i, j, t);
      // skip paired triangle
      if(t == pt) {
        continue;
      }
      mg[nVerts + nEdges + t].emplace_back(i);
    }
  }
  for(int i = 0; i < nTri; ++i) {
    const auto pT = this->getPairedCell(Cell{2, i}, triangulation);
    if(pT != -1) {
      mg[nVerts + nEdges + i].emplace_back(nVerts + nEdges + nTri + pT);
    }
    if(dim < 3) {
      continue;
    }
    // every other neighbor tetra is higher that the current triangle
    const auto nT = triangulation.getTriangleStarNumber(i);
    for(SimplexId j = 0; j < nT; ++j) {
      SimplexId T{};
      triangulation.getTriangleStar(i, j, T);
      // skip paired tetra
      if(T == pT) {
        continue;
      }
      mg[nVerts + nEdges + nTri + T].emplace_back(i);
    }
  }

  // fetch critical points
  std::vector<Cell> criticalPoints{};
  this->getCriticalPoints(criticalPoints, triangulation);
  std::sort(
    criticalPoints.begin(), criticalPoints.end(),
    [this, offset, &triangulation](const Cell &a, const Cell &b) {
      const auto oa = offset[this->getCellGreaterVertex(a, triangulation)];
      const auto ob = offset[this->getCellGreaterVertex(b, triangulation)];
      return oa < ob;
    });

  const auto cellId = [&](const Cell &c) {
    if(c.dim_ == 0)
      return c.id_;
    if(c.dim_ == 1)
      return c.id_ + nVerts;
    if(c.dim_ == 2)
      return c.id_ + nVerts + nEdges;
    if(c.dim_ == 3)
      return c.id_ + nVerts + nEdges + nTri;
    return -1;
  };

  // compare critical cells two by two
  for(size_t i = 0; i < criticalPoints.size() - 1; ++i) {
    mg[cellId(criticalPoints[i + 1])].emplace_back(cellId(criticalPoints[i]));
  }

  const auto sortedSimplices = this->topologicalSort(mg);
  if(sortedSimplices.empty()) {
    return {};
  }

  res.resize(sortedSimplices.size());
  for(size_t i = 0; i < sortedSimplices.size(); ++i) {
    res[sortedSimplices[i]] = sortedSimplices.size() - 1 - i;
  }

  this->printMsg("Computed Morse Discrete Function", 1.0, tm.getElapsedTime(),
                 this->threadNumber_);

  // check that discrete gradient is upheld
  for(int i = 0; i < nVerts; ++i) {
    const auto pe = this->getPairedCell(Cell{0, i}, triangulation);
    if(pe != -1 && res[i] < res[nVerts + pe]) {
      this->printErr("Vertex " + std::to_string(i)
                     + " should be higher than edge " + std::to_string(pe));
    }
  }
  for(int i = 0; i < nEdges; ++i) {
    const auto pt = this->getPairedCell(Cell{1, i}, triangulation);
    if(pt != -1 && res[nVerts + i] < res[nVerts + nEdges + pt]) {
      this->printErr("Edge " + std::to_string(i)
                     + " should be higher than triangle " + std::to_string(pt));
    }
  }
  for(int i = 0; i < nTri; ++i) {
    const auto pT = this->getPairedCell(Cell{2, i}, triangulation);
    if(pT != -1
       && res[nVerts + nEdges + i] < res[nVerts + nEdges + nTri + pT]) {
      this->printErr("Triangle " + std::to_string(i)
                     + " should be higher than tetra " + std::to_string(pT));
    }
  }

  return res;
}

template <typename triangulationType>
std::vector<ttk::DiscreteMorseFunction::Simplex>
  ttk::DiscreteMorseFunction::computeFiltrationOrder(
    const SimplexId *const offset,
    const triangulationType &triangulation) const {

  Timer tm{};

  const auto dim = triangulation.getDimensionality();

  const auto nVerts = triangulation.getNumberOfVertices();
  const auto nEdges = triangulation.getNumberOfEdges();
  const auto nTri = dim > 1 ? triangulation.getNumberOfTriangles() : 0;
  const auto nTetra = dim > 2 ? triangulation.getNumberOfCells() : 0;
  const auto num_simplices = nVerts + nEdges + nTri + nTetra;

  std::vector<Simplex> res(num_simplices);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    res[i].fillVert(i, offset);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nEdges; ++i) {
    const auto o = nVerts + i;
    res[o].fillEdge(i, o, offset, triangulation);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nTri; ++i) {
    const auto o = nVerts + nEdges + i;
    res[o].fillTriangle(i, o, offset, triangulation);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nTetra; ++i) {
    const auto o = nVerts + nEdges + nTri + i;
    res[o].fillTetra(i, o, offset, triangulation);
  }

  TTK_PSORT(this->threadNumber_, res.begin(), res.end());

  this->printMsg(
    "Computed filtration order", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return res;
}
