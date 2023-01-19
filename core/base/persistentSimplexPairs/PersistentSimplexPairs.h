/// \ingroup base
/// \class ttk::PersistentSimplexPairs
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date November 2021.
///
/// \brief Textbook algorithm to find persistence pairs.
///
/// This algorithm is described in "Algorithm and Theory of
/// Computation Handbook (Second Edition) - Special Topics and
/// Techniques" by Atallah and Blanton on page 97.

#pragma once

#include <AbstractTriangulation.h>
#include <Debug.h>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace ttk {

  class PersistentSimplexPairs : virtual public Debug {
  public:
    PersistentSimplexPairs();

    struct PersistencePair {
      /** first (lower) vertex id */
      SimplexId birth;
      /** second (higher) vertex id */
      SimplexId death;
      /** pair type (min-saddle: 0, saddle-saddle: 1, saddle-max: 2) */
      SimplexId type;

      PersistencePair(SimplexId b, SimplexId d, SimplexId t)
        : birth{b}, death{d}, type{t} {
      }
    };

    /**
     * @brief Preprocess all the required connectivity requests on the
     * triangulation.
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      if(data != nullptr) {
        const auto dim = data->getDimensionality();
        data->preconditionEdges();
        if(dim == 2) {
          data->preconditionCellEdges();
        } else if(dim == 3) {
          data->preconditionTriangles();
          data->preconditionTriangleEdges();
          data->preconditionCellTriangles();
        }
        this->nVerts_ = data->getNumberOfVertices();
        this->nEdges_ = data->getNumberOfEdges();
        this->nTri_ = dim > 1 ? data->getNumberOfTriangles() : 0;
        this->nTetra_ = dim > 2 ? data->getNumberOfCells() : 0;
      }
    }

    /**
     * @brief Compute the persistence pairs from the triangulation
     * simplicial complex
     */
    template <typename triangulationType>
    int computePersistencePairs(std::vector<PersistencePair> &pairs,
                                const SimplexId *const orderField,
                                const triangulationType &triangulation) const;

  private:
    /**
     * @brief Ad-hoc struct for sorting simplices
     */
    template <size_t n>
    struct Simplex {
      /** Index in the triangulation */
      SimplexId id_{};
      /** Order field value of the simplex vertices, sorted in
          decreasing order */
      std::array<SimplexId, n> vertsOrder_{};
      /** To compare two vertices according to the filtration (lexicographic
       * order) */
      friend bool operator<(const Simplex<n> &lhs, const Simplex<n> &rhs) {
        return lhs.vertsOrder_ < rhs.vertsOrder_;
      }
    };

    struct VertexSimplex : Simplex<1> {
      void fillVert(const SimplexId id, const SimplexId *const offsets) {
        this->id_ = id;
        this->vertsOrder_[0] = offsets[this->id_];
      }
    };
    struct EdgeSimplex : Simplex<2> {
      template <typename triangulationType>
      void fillEdge(const SimplexId id,
                    const SimplexId *const offsets,
                    const triangulationType &triangulation) {
        this->id_ = id;
        triangulation.getEdgeVertex(id, 0, this->vertsOrder_[0]);
        triangulation.getEdgeVertex(id, 1, this->vertsOrder_[1]);
        this->vertsOrder_[0] = offsets[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offsets[this->vertsOrder_[1]];
        // sort vertices in decreasing order
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }
    };
    struct TriangleSimplex : Simplex<3> {
      template <typename triangulationType>
      void fillTriangle(const SimplexId id,
                        const SimplexId *const offsets,
                        const triangulationType &triangulation) {
        this->id_ = id;
        triangulation.getTriangleVertex(id, 0, this->vertsOrder_[0]);
        triangulation.getTriangleVertex(id, 1, this->vertsOrder_[1]);
        triangulation.getTriangleVertex(id, 2, this->vertsOrder_[2]);
        this->vertsOrder_[0] = offsets[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offsets[this->vertsOrder_[1]];
        this->vertsOrder_[2] = offsets[this->vertsOrder_[2]];
        // sort vertices in decreasing order
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }
    };
    struct TetraSimplex : Simplex<4> {
      template <typename triangulationType>
      void fillTetra(const SimplexId id,
                     const SimplexId *const offsets,
                     const triangulationType &triangulation) {
        this->id_ = id;
        triangulation.getCellVertex(id, 0, this->vertsOrder_[0]);
        triangulation.getCellVertex(id, 1, this->vertsOrder_[1]);
        triangulation.getCellVertex(id, 2, this->vertsOrder_[2]);
        triangulation.getCellVertex(id, 3, this->vertsOrder_[3]);
        this->vertsOrder_[0] = offsets[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offsets[this->vertsOrder_[1]];
        this->vertsOrder_[2] = offsets[this->vertsOrder_[2]];
        this->vertsOrder_[3] = offsets[this->vertsOrder_[3]];
        // sort vertices in decreasing order
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }
    };

    template <typename triangulationType>
    void computeCellsOrder(std::vector<VertexSimplex> &verts,
                           std::vector<EdgeSimplex> &edges,
                           std::vector<TriangleSimplex> &triangles,
                           std::vector<TetraSimplex> &tetras,
                           std::array<std::vector<SimplexId>, 4> &cellsOrder,
                           const SimplexId *const offset,
                           const triangulationType &triangulation) const;

    inline SimplexId getCellId(const SimplexId cdim,
                               const SimplexId cid) const {
      if(cdim == 0) {
        return cid;
      } else if(cdim == 1) {
        return cid + this->nVerts_;
      } else if(cdim == 2) {
        return cid + this->nVerts_ + this->nEdges_;
      } else if(cdim == 3) {
        return cid + this->nVerts_ + this->nEdges_ + this->nTri_;
      }
      return -1;
    }

    template <typename triangulationType>
    int pairCells(std::vector<PersistencePair> &pairs,
                  std::vector<VertexSimplex> &verts,
                  std::vector<EdgeSimplex> &edges,
                  std::vector<TriangleSimplex> &triangles,
                  std::vector<TetraSimplex> &tetras,
                  const std::array<std::vector<SimplexId>, 4> &cellsOrder,
                  const triangulationType &triangulation) const;

    template <typename triangulationType,
              typename Container0,
              typename Container1>
    void
      pairCellsPerDim(std::vector<PersistencePair> &pairs,
                      const std::vector<Container0> &sortedFaces,
                      const std::vector<Container1> &sortedCells,
                      std::vector<bool> &isVisited,
                      std::vector<SimplexId> &partners,
                      std::array<std::vector<bool>, 4> &pairedSimplices,
                      const int dim,
                      const std::array<std::vector<SimplexId>, 4> &cellsOrder,
                      const triangulationType &triangulation) const {

      isVisited.resize(sortedFaces.size(), false);
      partners.resize(sortedFaces.size());
      std::fill(partners.begin(), partners.end(), -1);
      const auto cmp
        = [&cellsOrder, &dim](const SimplexId a, const SimplexId b) {
            return cellsOrder[dim - 1][a] > cellsOrder[dim - 1][b];
          };
      using Container = std::set<SimplexId, decltype(cmp)>;
      std::vector<Container> boundaries(sortedCells.size(), Container(cmp));

      for(size_t j = 0; j < sortedCells.size(); ++j) {
        const auto &c{sortedCells[j]};
        const auto tau = this->eliminateBoundaries(
          c.id_, dim, isVisited, boundaries, partners, triangulation);
        if(tau != -1) {
          const auto &pc{sortedFaces[cellsOrder[dim - 1][tau]]};
          pairedSimplices[dim - 1][pc.id_] = true;
          pairedSimplices[dim][c.id_] = true;
        }
      }

      for(size_t i = 0; i < partners.size(); ++i) {
        if(partners[i] == -1) {
          continue;
        }
        const auto &c{sortedFaces[cellsOrder[dim - 1][i]]};
        const auto &pc{sortedCells[cellsOrder[dim][partners[i]]]};
        // only record pairs with non-null persistence
        if(c.vertsOrder_[0] != pc.vertsOrder_[0]) {
          pairs.emplace_back(i, partners[i], dim - 1);
        }
      }
    }

    template <typename triangulationType, typename Container>
    SimplexId
      eliminateBoundaries(const SimplexId c,
                          const int dim,
                          std::vector<bool> &isVisited,
                          std::vector<Container> &boundaries,
                          std::vector<SimplexId> &partners,
                          const triangulationType &triangulation) const {

      auto &boundary{boundaries[c]};
      const auto expandBoundary = [&boundary, &isVisited](const SimplexId s) {
        if(!isVisited[s]) {
          boundary.emplace(s);
          isVisited[s] = true;
        } else {
          const auto it{boundary.find(s)};
          boundary.erase(it);
          isVisited[s] = false;
        }
      };

      const auto clearBoundary = [&boundary, &isVisited]() {
        for(const auto e : boundary) {
          isVisited[e] = false;
        }
      };

      const auto addBoundaryEl = [&triangulation, &dim, &expandBoundary](
                                   const SimplexId a, const int lid) {
        SimplexId s{};
        if(dim == 1) {
          triangulation.getEdgeVertex(a, lid, s);
        } else if(dim == 2) {
          triangulation.getTriangleEdge(a, lid, s);
        } else if(dim == 3) {
          triangulation.getCellTriangle(a, lid, s);
        }
        expandBoundary(s);
      };

      const auto addBoundary = [&addBoundaryEl, &dim](const SimplexId a) {
        for(int i = 0; i < dim + 1; ++i) {
          addBoundaryEl(a, i);
        }
      };
      addBoundary(c);

      while(!boundary.empty()) {
        // youngest cell on boundary
        const auto tau{*boundary.begin()};
        const auto partnerTau{partners[tau]};
        if(partnerTau == -1) {
          partners[tau] = c;
          clearBoundary();
          return tau;
        }
        if(boundaries[partnerTau].empty()) {
          addBoundary(partnerTau);
        } else {
          // merge boundaries
          for(const auto s : boundaries[partnerTau]) {
            expandBoundary(s);
          }
        }
      }

      clearBoundary();
      return -1;
    }

    SimplexId nVerts_{0};
    SimplexId nEdges_{0};
    SimplexId nTri_{0};
    SimplexId nTetra_{0};
  };

} // namespace ttk

template <typename triangulationType>
int ttk::PersistentSimplexPairs::pairCells(
  std::vector<PersistencePair> &pairs,
  std::vector<VertexSimplex> &verts,
  std::vector<EdgeSimplex> &edges,
  std::vector<TriangleSimplex> &triangles,
  std::vector<TetraSimplex> &tetras,
  const std::array<std::vector<SimplexId>, 4> &cellsOrder,
  const triangulationType &triangulation) const {

  const auto dim{triangulation.getDimensionality()};
  std::array<std::vector<PersistencePair>, 3> pairsPerDim{};

  std::array<std::vector<bool>, 4> pairedSimplices{};
  std::array<std::vector<bool>, 3> isVisited{};
  pairedSimplices[0].resize(this->nVerts_, false);
  isVisited[0].resize(this->nVerts_, false);
  if(dim > 0) {
    pairedSimplices[1].resize(this->nEdges_, false);
    isVisited[1].resize(this->nEdges_, false);
  }
  if(dim > 1) {
    pairedSimplices[2].resize(this->nTri_, false);
    isVisited[2].resize(this->nTri_, false);
  }
  if(dim > 2) {
    pairedSimplices[3].resize(this->nTetra_, false);
  }

#pragma omp parallel num_threads(this->threadNumber_)
#pragma omp sections
  {
#pragma omp section
    {
      Timer tm{};
      std::vector<SimplexId> partners{};
      this->pairCellsPerDim(pairsPerDim[0], verts, edges, isVisited[0],
                            partners, pairedSimplices, 1, cellsOrder,
                            triangulation);
      this->printMsg("Computed " + std::to_string(pairsPerDim[0].size())
                       + " pairs of dimension 0",
                     1.0, tm.getElapsedTime(), 1);
    }
#pragma omp section
    if(dim > 1) {
      Timer tm{};
      std::vector<SimplexId> partners{};
      this->pairCellsPerDim(pairsPerDim[1], edges, triangles, isVisited[1],
                            partners, pairedSimplices, 2, cellsOrder,
                            triangulation);
      this->printMsg("Computed " + std::to_string(pairsPerDim[1].size())
                       + " pairs of dimension 1",
                     1.0, tm.getElapsedTime(), 1);
    }
#pragma omp section
    if(dim > 2) {
      Timer tm{};
      std::vector<SimplexId> partners{};
      this->pairCellsPerDim(pairsPerDim[2], triangles, tetras, isVisited[2],
                            partners, pairedSimplices, 3, cellsOrder,
                            triangulation);
      this->printMsg("Computed " + std::to_string(pairsPerDim[2].size())
                       + " pairs of dimension 2",
                     1.0, tm.getElapsedTime(), 1);
    }
  }

  pairs = std::move(pairsPerDim[0]);
  pairs.insert(pairs.end(), pairsPerDim[2].begin(), pairsPerDim[2].end());
  pairs.insert(pairs.end(), pairsPerDim[1].begin(), pairsPerDim[1].end());

  Timer tm{};
  const auto nPairs{pairs.size()};
  for(size_t i = 0; i < pairedSimplices.size(); ++i) {
    for(size_t j = 0; j < pairedSimplices[i].size(); ++j) {
      if(!pairedSimplices[i][j]) {
        pairs.emplace_back(j, -1, i);
      }
    }
  }
  this->printMsg(
    "Detected " + std::to_string(pairs.size() - nPairs) + " infinite pairs",
    1.0, tm.getElapsedTime(), 1);

  return 0;
}

template <typename triangulationType>
int ttk::PersistentSimplexPairs::computePersistencePairs(
  std::vector<ttk::PersistentSimplexPairs::PersistencePair> &pairs,
  const SimplexId *const orderField,
  const triangulationType &triangulation) const {

  Timer tm{};

  std::array<std::vector<SimplexId>, 4> cellsOrder{};
  std::vector<VertexSimplex> verts{};
  std::vector<EdgeSimplex> edges{};
  std::vector<TriangleSimplex> triangles{};
  std::vector<TetraSimplex> tetras{};

  this->computeCellsOrder(
    verts, edges, triangles, tetras, cellsOrder, orderField, triangulation);
  this->pairCells(
    pairs, verts, edges, triangles, tetras, cellsOrder, triangulation);

  this->printMsg("Computed " + std::to_string(pairs.size())
                   + " persistence pair" + (pairs.size() > 1 ? "s" : ""),
                 1.0, tm.getElapsedTime(), 1);

  return 0;
}

template <typename triangulationType>
void ttk::PersistentSimplexPairs::computeCellsOrder(
  std::vector<VertexSimplex> &verts,
  std::vector<EdgeSimplex> &edges,
  std::vector<TriangleSimplex> &triangles,
  std::vector<TetraSimplex> &tetras,
  std::array<std::vector<SimplexId>, 4> &cellsOrder,
  const SimplexId *const offsets,
  const triangulationType &triangulation) const {

  Timer tm{};

  cellsOrder[0].resize(this->nVerts_);
  cellsOrder[1].resize(this->nEdges_);
  cellsOrder[2].resize(this->nTri_);
  cellsOrder[3].resize(this->nTetra_);
  verts.resize(this->nVerts_);
  edges.resize(this->nEdges_);
  triangles.resize(this->nTri_);
  tetras.resize(this->nTetra_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->nVerts_; ++i) {
      verts[i].fillVert(i, offsets);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->nEdges_; ++i) {
      edges[i].fillEdge(i, offsets, triangulation);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->nTri_; ++i) {
      triangles[i].fillTriangle(i, offsets, triangulation);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->nTetra_; ++i) {
      tetras[i].fillTetra(i, offsets, triangulation);
    }
  }

  TTK_PSORT(this->threadNumber_, verts.begin(), verts.end());
  TTK_PSORT(this->threadNumber_, edges.begin(), edges.end());
  TTK_PSORT(this->threadNumber_, triangles.begin(), triangles.end());
  TTK_PSORT(this->threadNumber_, tetras.begin(), tetras.end());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < verts.size(); ++i) {
      cellsOrder[0][verts[i].id_] = i;
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < edges.size(); ++i) {
      cellsOrder[1][edges[i].id_] = i;
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < triangles.size(); ++i) {
      cellsOrder[2][triangles[i].id_] = i;
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < tetras.size(); ++i) {
      cellsOrder[3][tetras[i].id_] = i;
    }
  }

  this->printMsg(
    "Computed filtration order", 1.0, tm.getElapsedTime(), this->threadNumber_);
}
