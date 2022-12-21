#pragma once

#include <CinemaImagingNative.h>

class vtkMultiBlockDataSet;
class vtkPointSet;

class ttkCinemaImagingNative : public ttk::CinemaImagingNative {
public:
  ttkCinemaImagingNative();
  ~ttkCinemaImagingNative() override;

  int RenderVTKObject(vtkMultiBlockDataSet *outputImages,

                      vtkPointSet *inputObject,
                      vtkPointSet *inputGrid) const;
};
