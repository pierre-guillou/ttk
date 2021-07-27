/// \ingroup vtk
/// \class ttkDiscreteMorseFunction
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date August 2021
///
/// \brief TTK VTK-filter that computes the Discrete Morse Function.
///
/// This VTK filter uses the ttk::dcg::DiscreteMorseFunction module to compute
/// the Discrete Morse Function of an input scalar field defined on
/// the vertices of a Simplicial Complex.
///
/// \param Input vtkUnstructuredGrid.
/// \param Output vtkUnstructuredGrid.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::dcg::DiscreteMorseFunction
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkDiscreteMorseFunctionModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <DiscreteMorseFunction.h>

class TTKDISCRETEMORSEFUNCTION_EXPORT ttkDiscreteMorseFunction
  : public ttkAlgorithm,
    protected ttk::DiscreteMorseFunction {

public:
  vtkSetMacro(DMFArrayName, const std::string &);
  vtkGetMacro(DMFArrayName, std::string);

  vtkSetMacro(FiltrArrayName, const std::string &);
  vtkGetMacro(FiltrArrayName, std::string);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  static ttkDiscreteMorseFunction *New();
  vtkTypeMacro(ttkDiscreteMorseFunction, ttkAlgorithm);

protected:
  ttkDiscreteMorseFunction();
  ~ttkDiscreteMorseFunction() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::vector<SimplexId> dmf{}, filtr{};
  std::string DMFArrayName{"DiscreteMorseFunction"};
  std::string FiltrArrayName{"FiltrationOrder"};
  bool ForceInputOffsetScalarField{false};
};
