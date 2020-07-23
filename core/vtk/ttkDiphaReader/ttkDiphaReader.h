/// \ingroup vtk
/// \class ttkDiphaReader
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date January 2021.
/// \brief ttkDiphaReader - Object File Format Reader
///
/// Load a .dipha file into VTK format

#pragma once

#include <Debug.h>
#include <ttkDiphaReaderModule.h>

#include <vtkUnstructuredGridAlgorithm.h>

class TTKDIPHAREADER_EXPORT ttkDiphaReader
  : public vtkUnstructuredGridAlgorithm,
    protected ttk::Debug {
public:
  vtkTypeMacro(ttkDiphaReader, vtkUnstructuredGridAlgorithm);

  static ttkDiphaReader *New();

  void PrintSelf(std::ostream &os, vtkIndent indent) override;

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  ttkDiphaReader();

  int readPersistenceDiagram(std::ifstream &, vtkUnstructuredGrid *) const;
  int readImageData(std::ifstream &, vtkUnstructuredGrid *) const;

  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *) override;

private:
  ttkDiphaReader(const ttkDiphaReader &) = delete;
  void operator=(const ttkDiphaReader &) = delete;

  char *FileName{};
};
