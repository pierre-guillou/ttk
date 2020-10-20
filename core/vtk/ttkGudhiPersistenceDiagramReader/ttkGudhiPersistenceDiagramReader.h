/// \ingroup vtk
/// \class ttkGudhiPersistenceDiagramReader
/// \author Charles Gueunet <charles.gueunet@kitware.com>
/// \date December 2017.
/// \brief ttkGudhiPersistenceDiagramReader - Gudhi Persistence Diagram Reader
///
/// Load a .pers file into VTK format
///
/// Note: This reader is not able to deal with comment on the file

#pragma once

#include <Debug.h>
#include <ttkGudhiPersistenceDiagramReaderModule.h>

#include <vtkUnstructuredGridAlgorithm.h>

class TTKGUDHIPERSISTENCEDIAGRAMREADER_EXPORT ttkGudhiPersistenceDiagramReader
  : public vtkUnstructuredGridAlgorithm,
    protected ttk::Debug {
public:
  vtkTypeMacro(ttkGudhiPersistenceDiagramReader, vtkUnstructuredGridAlgorithm);

  static ttkGudhiPersistenceDiagramReader *New();

  void PrintSelf(std::ostream &os, vtkIndent indent) override;

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  ttkGudhiPersistenceDiagramReader();

  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *) override;

private:
  ttkGudhiPersistenceDiagramReader(const ttkGudhiPersistenceDiagramReader &)
    = delete;
  void operator=(const ttkGudhiPersistenceDiagramReader &) = delete;

  char *FileName{};
};
