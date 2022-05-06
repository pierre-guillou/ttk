/// \ingroup base
/// \class ttk::Mapper
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date May 2022.
///
/// \brief TTK processing package for scalar field smoothing.
///
/// This class is a dummy example for the development of TTK classes. It
/// smooths an input scalar field by averaging the scalar values on the link
/// of each vertex.
///
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkMapper.cpp %for a usage example.

#pragma once

#include <Triangulation.h>

namespace ttk {

  class Mapper : virtual public Debug {
  public:
    Mapper();

    inline void setDimensionNumber(const int &dimensionNumber) {
      dimensionNumber_ = dimensionNumber;
    }

    inline void setInputDataPointer(void *data) {
      inputData_ = data;
    }

    inline void setOutputDataPointer(void *data) {
      outputData_ = data;
    }

    inline void setMaskDataPointer(const char *const mask) {
      mask_ = mask;
    }

    inline void
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation != nullptr) {
        triangulation->preconditionVertexNeighbors();
      }
    }

    template <class dataType, class triangulationType>
    int execute(const triangulationType *triangulation,
                const int &numberOfIterations) const;

  protected:
    int dimensionNumber_{1};
    void *inputData_{nullptr}, *outputData_{nullptr};
    const char *mask_{nullptr};
  };

} // namespace ttk

// template functions
template <class dataType, class triangulationType>
int ttk::Mapper::execute(const triangulationType *triangulation,
                         const int &numberOfIterations) const {

  Timer tm{};

  printMsg("Done!", 1, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}
