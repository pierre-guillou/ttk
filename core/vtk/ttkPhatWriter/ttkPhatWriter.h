/// \ingroup vtk
/// \class ttkPhatWriter
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date July 2020
/// \brief ttkPhatWriter - Dipha Image Data Format Writer
///
/// Writes a Dipha Cubical Complex file from a VTK Image Data or a
/// Dipha Explicit Complex from a VTK Unstructured Grid dataset

#pragma once

#include <ttkAlgorithm.h>
#include <ttkPhatWriterModule.h>

#include <fstream>

class TTKPHATWRITER_EXPORT ttkPhatWriter : public ttkAlgorithm {

public:
  vtkTypeMacro(ttkPhatWriter, ttkAlgorithm);

  static ttkPhatWriter *New();

  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

  // expose vtkWriter methods (duck-typing)
  int Write();
  vtkDataObject *GetInput();
  void SetInputData(vtkDataObject *input);

protected:
  // Regular writer management.
  ttkPhatWriter();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int writeUnstructuredGrid(vtkDataObject *);

  int OpenFile();

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkPhatWriter(const ttkPhatWriter &) = delete;
  void operator=(const ttkPhatWriter &) = delete;
};
