#pragma once

#include <ttkAlgorithmModule.h>

#include <Debug.h>
#include <array>
#include <unordered_map>
#include <vtkType.h>

class vtkDataSet;
class vtkImageData;
class vtkPointSet;
class vtkPoints;
class vtkCellArray;
namespace ttk {
  class AbstractTriangulation;
}

// allow to store 2 triangulation pointers (for implicit + periodic)
using RegistryTriangulation
  = std::array<std::unique_ptr<ttk::AbstractTriangulation>, 2>;

struct RegistryValue {
  RegistryTriangulation triangulation;
  vtkDataSet *owner;

  vtkMTimeType cellModTime{0};

  int extent[6];
  double origin[3];
  double spacing[3];
  int dimensions[3];

  RegistryValue(vtkDataSet *dataSet,
                ttk::AbstractTriangulation *triangulation0_,
                ttk::AbstractTriangulation *triangulation1_ = nullptr);
  bool isValid(vtkDataSet *dataSet) const;
};

using RegistryKey = long long;
using Registry = std::unordered_map<RegistryKey, RegistryValue>;

class TTKALGORITHM_EXPORT ttkTriangulationFactory : public ttk::Debug {
public:
  static ttk::AbstractTriangulation *GetTriangulation(int debugLevel,
                                                      vtkDataSet *object);
  static int
    SwitchToPeriodicTriangulation(ttk::AbstractTriangulation *triangulation);

  static ttkTriangulationFactory Instance;
  static RegistryKey GetKey(vtkDataSet *dataSet);

#ifdef _WIN32
  // to fix a weird MSVC warning about unique_ptr inside
  // unordered_map, this dummy class member should be declared before
  // the Registry
  RegistryTriangulation dummy{};
#endif // _WIN32
  Registry registry;

private:
  RegistryTriangulation CreateImplicitTriangulation(vtkImageData *image);
  RegistryTriangulation CreateExplicitTriangulation(vtkPointSet *pointSet);
  RegistryTriangulation CreateTriangulation(vtkDataSet *dataSet);
  int FindImplicitTriangulation(ttk::AbstractTriangulation *&triangulation,
                                vtkImageData *image);

  ttkTriangulationFactory();
};
