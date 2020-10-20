/// \ingroup vtk
/// \class ttkPerseusWriter
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date July 2020
/// \brief ttkPerseusWriter - Perseus Cubical Grid Format Writer
///
/// Writes a Perseus Cubical Grid file from a VTK Image Data.

#pragma once

#include <vtkDataSetWriter.h>

#include <ttkAlgorithm.h>
#include <ttkPerseusWriterModule.h>

#include <fstream>

class TTKPERSEUSWRITER_EXPORT ttkPerseusWriter : public ttkAlgorithm {

public:
  vtkTypeMacro(ttkPerseusWriter, ttkAlgorithm);

  static ttkPerseusWriter *New();

  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

  // expose vtkWriter methods (duck-typing)
  int Write();
  vtkDataObject *GetInput();
  void SetInputData(vtkDataObject *input);

protected:
  // Regular writer management
  ttkPerseusWriter();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int writeImageData(vtkDataObject *);
  int writeUnstructuredGrid(vtkDataObject *);

  int OpenFile();

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkPerseusWriter(const ttkPerseusWriter &) = delete;
  void operator=(const ttkPerseusWriter &) = delete;
};
