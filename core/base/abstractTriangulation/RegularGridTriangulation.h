/// \ingroup base
/// \class ttk::RegularGridTriangulation
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date December 2022.
///
/// \brief RegularGridTriangulation is an abstract subclass of
/// ttk::AbstractTriangulation that exposes a common API for
/// triangulation on regular grids. This class is meant to be
/// implemented in ttk::ImplicitTriangulation and
/// ttk::PeriodicImplicitTriangulation.
///
/// \sa ttk::Triangulation

#pragma once

#include <AbstractTriangulation.h>

namespace ttk {

  class ImplicitTriangulation;
  class PeriodicImplicitTriangulation;

  class RegularGridTriangulation : public AbstractTriangulation {
    friend class ttk::ImplicitTriangulation;
    friend class ttk::PeriodicImplicitTriangulation;

  public:
    ~RegularGridTriangulation() override = default;

    virtual int setInputGrid(const float &xOrigin,
                             const float &yOrigin,
                             const float &zOrigin,
                             const float &xSpacing,
                             const float &ySpacing,
                             const float &zSpacing,
                             const SimplexId &xDim,
                             const SimplexId &yDim,
                             const SimplexId &zDim)
      = 0;

  protected:
    virtual void vertexToPosition2d(const SimplexId vertex,
                                    SimplexId p[2]) const = 0;
    virtual void vertexToPosition(const SimplexId vertex,
                                  SimplexId p[3]) const = 0;
    virtual void triangleToPosition2d(const SimplexId triangle,
                                      SimplexId p[2]) const = 0;
    virtual void triangleToPosition(const SimplexId triangle,
                                    const int k,
                                    SimplexId p[3]) const = 0;
    virtual void tetrahedronToPosition(const SimplexId tetrahedron,
                                       SimplexId p[3]) const = 0;

  private:
    SimplexId findEdgeFromVertices(const SimplexId v0,
                                   const SimplexId v1) const;
    SimplexId findTriangleFromVertices(std::array<SimplexId, 3> &verts) const;
  };

} // namespace ttk
