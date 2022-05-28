#include <SurfaceGeometrySmoother.h>

#include <array>
#include <iostream>

namespace ttk {

  class Quadrangulation : virtual public Debug {
  public:
    int preconditionVertexNeighbors();

    inline void setInputCells(const SimplexId cellNumber,
                              void *const quadCells) {
      this->nCells_ = cellNumber;
      this->cells_ = static_cast<const Quad *>(quadCells);
    }

    inline void setInputPoints(const SimplexId pointNumber,
                               void *const pointCoords) {
      this->nVerts_ = pointNumber;
      this->vertCoords_ = static_cast<const Point *>(pointCoords);
    }

    inline int isVertexExtraordinary(const SimplexId v) {
      const size_t ordinary_valence{4};
      return this->vertexNeighbors_[v].size() != ordinary_valence;
    }

    inline void
      getVertexPoint(const SimplexId v, float &x, float &y, float &z) const {
      const auto &c{this->vertCoords_[v]};
      x = c[0];
      y = c[1];
      z = c[2];
    }

    inline void
      getCellVertex(const SimplexId c, const int l, SimplexId &v) const {
      v = this->cells_[c][l];
    }

    inline SimplexId getVertexNeighborNumber(const SimplexId v) const {
      return this->vertexNeighbors_[v].size();
    }
    inline void
      getVertexNeighbor(const SimplexId v, const int l, SimplexId &n) const {
      n = this->vertexNeighbors_[v][l];
    }

    inline SimplexId getCellVertexNumber(const SimplexId ttkNotUsed(c)) const {
      return 4;
    }
    inline int getDimensionality() const {
      return 2;
    }
    inline SimplexId getNumberOfVertices() const {
      return this->nVerts_;
    }
    inline SimplexId getNumberOfCells() const {
      return this->nCells_;
    }

  private:
    /**
     * @brief Ad-hoc quad data structure (4 vertex ids)
     */
    using Quad = std::array<ttk::LongSimplexId, 4>;
    /**
     * @brief Ad-hoc vertex coordinates data structure (3 floats)
     */
    using Point = ttk::SurfaceGeometrySmoother::Point;

    const Point *vertCoords_{};
    const Quad *cells_{};
    SimplexId nVerts_{};
    SimplexId nCells_{};
    FlatJaggedArray vertexNeighbors_{};
  };

} // namespace ttk
