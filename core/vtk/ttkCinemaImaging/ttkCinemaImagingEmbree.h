#pragma once

#include <CinemaImagingEmbree.h>

class vtkMultiBlockDataSet;
class vtkPointSet;

class ttkCinemaImagingEmbree : public ttk::CinemaImagingEmbree {
public:
  ttkCinemaImagingEmbree();
  ~ttkCinemaImagingEmbree() override;

  int RenderVTKObject(vtkMultiBlockDataSet *outputImages,

                      vtkPointSet *inputObject,
                      vtkPointSet *inputGrid) const;
};
