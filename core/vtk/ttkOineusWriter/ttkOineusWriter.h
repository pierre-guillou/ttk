/// \ingroup vtk
/// \class ttkOineusWriter
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date October 2022
/// \brief ttkOineusWriter - Oineus Format Writer

#pragma once

#include <ttkAlgorithm.h>
#include <ttkOineusWriterModule.h>

#include <fstream>

class TTKOINEUSWRITER_EXPORT ttkOineusWriter : public ttkAlgorithm {

public:
  vtkTypeMacro(ttkOineusWriter, ttkAlgorithm);

  static ttkOineusWriter *New();

  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

  // expose vtkWriter methods (duck-typing)
  int Write();
  vtkDataObject *GetInput();
  void SetInputData(vtkDataObject *input);

protected:
  // Regular writer management.
  ttkOineusWriter();
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int writeUnstructuredGrid(vtkDataObject *);

  int OpenFile();

  char *Filename{};
  std::ofstream Stream{};

private:
  ttkOineusWriter(const ttkOineusWriter &) = delete;
  void operator=(const ttkOineusWriter &) = delete;
};
