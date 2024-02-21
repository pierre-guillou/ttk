/// \ingroup vtk
/// \class ttkAbsDiff
/// \author Pierre Guillou <pierre.guillou@mines-paris.org>
/// \date 21/02/2024
///
/// \brief TTK VTK-filter that compute the absolute difference between
/// all scalar fields of two datasets

#pragma once

// VTK Module
#include <ttkAbsDiffModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKABSDIFF_EXPORT ttkAbsDiff : public ttkAlgorithm {
public:
  static ttkAbsDiff *New();
  vtkTypeMacro(ttkAbsDiff, ttkAlgorithm);

protected:
  ttkAbsDiff();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
