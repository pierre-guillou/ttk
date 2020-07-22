/// \ingroup vtk
/// \class ttkDiphaWriter
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date July 2020
/// \brief ttkDiphaWriter - Dipha Image Data Format Writer
///
/// Writes a Dipha Cubical Complex file from a VTK Image Data or a
/// Dipha Explicit Complex from a VTK Unstructured Grid dataset

#pragma once

#include <ttkAlgorithm.h>
#include <ttkDiphaWriterModule.h>

#include <fstream>

class TTKDIPHAWRITER_EXPORT ttkDiphaWriter : public ttkAlgorithm {

public:
  vtkTypeMacro(ttkDiphaWriter, ttkAlgorithm);

  static ttkDiphaWriter *New();

  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

  // expose vtkWriter methods (duck-typing)
  int Write();
  vtkDataObject *GetInput();
  void SetInputData(vtkDataObject *input);

protected:
  // Regular writer management.
  ttkDiphaWriter();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int writeImageData(vtkDataObject *);
  int writeUnstructuredGrid(vtkDataObject *);

  int OpenFile();

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkDiphaWriter(const ttkDiphaWriter &) = delete;
  void operator=(const ttkDiphaWriter &) = delete;
};
