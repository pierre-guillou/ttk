/// \ingroup vtk
/// \class ttkEireneWriter
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date July 2020
/// \brief ttkEireneWriter - Dipha Image Data Format Writer
///
/// Writes a Dipha Cubical Complex file from a VTK Image Data or a
/// Dipha Explicit Complex from a VTK Unstructured Grid dataset

#pragma once

#include <ttkAlgorithm.h>
#include <ttkEireneWriterModule.h>

#include <fstream>

class TTKEIRENEWRITER_EXPORT ttkEireneWriter : public ttkAlgorithm {

public:
  vtkTypeMacro(ttkEireneWriter, ttkAlgorithm);

  static ttkEireneWriter *New();

  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

  // expose vtkWriter methods (duck-typing)
  int Write();
  vtkDataObject *GetInput();
  void SetInputData(vtkDataObject *input);

protected:
  // Regular writer management.
  ttkEireneWriter();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int writeUnstructuredGrid(vtkDataObject *);

  int OpenFile();

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkEireneWriter(const ttkEireneWriter &) = delete;
  void operator=(const ttkEireneWriter &) = delete;
};
