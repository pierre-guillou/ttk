/// \ingroup vtk
/// \class ttkSimplicialComplexWriter
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date July 2020
/// \brief ttkSimplicialComplexWriter - Dipha Image Data Format Writer
///
/// Writes a Dipha Cubical Complex file from a VTK Image Data or a
/// Dipha Explicit Complex from a VTK Unstructured Grid dataset

#pragma once

#include <ttkAlgorithm.h>
#include <ttkSimplicialComplexWriterModule.h>

#include <fstream>

class TTKSIMPLICIALCOMPLEXWRITER_EXPORT ttkSimplicialComplexWriter
  : public ttkAlgorithm {

public:
  vtkTypeMacro(ttkSimplicialComplexWriter, ttkAlgorithm);

  static ttkSimplicialComplexWriter *New();

  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

  // expose vtkWriter methods (duck-typing)
  int Write();
  vtkDataObject *GetInput();
  void SetInputData(vtkDataObject *input);

protected:
  // Regular writer management.
  ttkSimplicialComplexWriter();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int writeUnstructuredGrid(vtkDataObject *);

  int OpenFile();

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkSimplicialComplexWriter(const ttkSimplicialComplexWriter &) = delete;
  void operator=(const ttkSimplicialComplexWriter &) = delete;
};
