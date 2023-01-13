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
#include <VisitedMask.h>

#include <algorithm>
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
    void computeCellsOrder(std::array<std::vector<SimplexId>, 4> &cellsOrder,
                           const SimplexId *const offset,
                           const triangulationType &triangulation) const;

    inline void addCellBoundary(const Simplex &c, VisitedMask &boundary) const {
      for(SimplexId i = 0; i < c.dim_ + 1; ++i) {
        const auto f{c.faceIds_[i]};
        if(!boundary.isVisited_[f]) {
          boundary.insert(f);
        } else {
          boundary.remove(f);
        }
      }
    }

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

    SimplexId
      eliminateBoundaries(const Simplex &c,
                          VisitedMask &boundary,
                          std::vector<SimplexId> &partners,
                          const std::vector<Simplex> &filtration,
                          const std::vector<SimplexId> &filtOrder) const;

    int pairCells(std::vector<PersistencePair> &pairs,
                  const std::vector<Simplex> &filtration,
                  const std::vector<SimplexId> &filtOrder) const;

    SimplexId nVerts_{0};
    SimplexId nEdges_{0};
    SimplexId nTri_{0};
    SimplexId nTetra_{0};
  };

} // namespace ttk

template <typename triangulationType>
int ttk::PersistentSimplexPairs::computePersistencePairs(
  std::vector<ttk::PersistentSimplexPairs::PersistencePair> &pairs,
  const SimplexId *const orderField,
  const triangulationType &triangulation) const {

  Timer tm{};

  // every simplex in the triangulation, sorted by filtration
  const auto filtration
    = this->computeFiltrationOrder(orderField, triangulation);

  // simplex id -> filtration order
  std::vector<SimplexId> filtOrder(filtration.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < filtration.size(); ++i) {
    filtOrder[filtration[i].cellId_] = i;
  }

  this->pairCells(pairs, filtration, filtOrder);

  this->printMsg("Computed " + std::to_string(pairs.size())
                   + " persistence pair" + (pairs.size() > 1 ? "s" : ""),
                 1.0, tm.getElapsedTime(), 1);

  return 0;
}

template <typename triangulationType>
void ttk::PersistentSimplexPairs::computeCellsOrder(
  std::array<std::vector<SimplexId>, 4> &cellsOrder,
  const SimplexId *const offsets,
  const triangulationType &triangulation) const {

  Timer tm{};

  cellsOrder[0].resize(this->nVerts_);
  cellsOrder[1].resize(this->nEdges_);
  cellsOrder[2].resize(this->nTri_);
  cellsOrder[3].resize(this->nTetra_);

  std::vector<VertexSimplex> verts(this->nVerts_);
  std::vector<EdgeSimplex> edges(this->nEdges_);
  std::vector<TriangleSimplex> triangles(this->nTri_);
  std::vector<TetraSimplex> tetras(this->nTetra_);

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
