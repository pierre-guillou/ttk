/// \ingroup vtk
/// \class ttkMapper
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date May 2022.
///
/// \brief TTK Mapper VTK-filter.
///
/// This class generates a mapper from a data-set with a given number
/// of buckets. Data-set vertices and edges are placed into buckets
/// from their input scalar field value then connected components are
/// detected per bucket.
///
/// \param Input vtkDataSet
/// \param Output Nodes, Arcs (both vtkUnstructuredGrid), Segmentation
/// (vtkDataSet)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).

#pragma once

// VTK includes

// VTK Module
#include <ttkMapperModule.h>

// ttk code includes
#include <Mapper.h>
#include <ttkAlgorithm.h>

class TTKMAPPER_EXPORT ttkMapper : public ttkAlgorithm, protected ttk::Mapper {
public:
  static ttkMapper *New();

  vtkTypeMacro(ttkMapper, ttkAlgorithm);

  vtkSetMacro(NumberOfBuckets, int);
  vtkGetMacro(NumberOfBuckets, int);

  vtkSetMacro(ComputeHighDimBarycenters, bool);
  vtkGetMacro(ComputeHighDimBarycenters, bool);

  void SetHighDimCoords(const std::string &s) {
    HighDimCoords.push_back(s);
    Modified();
  }

  void ClearHighDimCoords() {
    HighDimCoords.clear();
    Modified();
  }

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

protected:
  ttkMapper();

  ~ttkMapper() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  std::vector<std::string> HighDimCoords{};
};
