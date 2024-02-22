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
#include <ttkMacros.h>

class TTKABSDIFF_EXPORT ttkAbsDiff : public ttkAlgorithm {
public:
  enum class FUNCTION {
    ABSDIFF = 0,
    RATIO = 1,
    RELDIFF = 2,
  };

  static ttkAbsDiff *New();
  vtkTypeMacro(ttkAbsDiff, ttkAlgorithm);

  ttkSetEnumMacro(Function, FUNCTION);
  vtkGetEnumMacro(Function, FUNCTION);

protected:
  ttkAbsDiff();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <typename ArrT, typename Func>
  void compute(ArrT res, const ArrT arr0, const ArrT arr1, const Func &f);

  template <typename ArrT>
  void dispatch(ArrT res, const ArrT arr0, const ArrT arr1);

  FUNCTION Function{FUNCTION::ABSDIFF};
};
