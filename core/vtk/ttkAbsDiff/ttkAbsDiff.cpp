#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include "ttkAbsDiff.h"
#include <ttkUtils.h>

vtkStandardNewMacro(ttkAbsDiff);

ttkAbsDiff::ttkAbsDiff() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkAbsDiff::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkAbsDiff::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkAbsDiff::RequestData(vtkInformation *ttkNotUsed(request),
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector) {

  this->setDebugLevel(2);

  const auto nInputs = inputVector[0]->GetNumberOfInformationObjects();
  if(nInputs < 2) {
    this->printErr("Need 2 inputs, got " + std::to_string(nInputs));
    return 0;
  }

  if(nInputs > 2) {
    this->printWrn("Got " + std::to_string(nInputs)
                   + " inputs, only the first two will be considered");
  }

  const auto input0 = vtkDataSet::GetData(inputVector[0], 0);
  const auto input1 = vtkDataSet::GetData(inputVector[0], 1);
  auto output = vtkDataSet::GetData(outputVector);

  if(input0->GetNumberOfPoints() != input1->GetNumberOfPoints()) {
    this->printErr("Input datasets should have the same number of points");
    return 0;
  }

  if(input0->GetNumberOfCells() != input1->GetNumberOfCells()) {
    this->printErr("Input datasets should have the same number of cells");
    return 0;
  }

  if(input0->GetPointData()->GetNumberOfArrays()
     != input1->GetPointData()->GetNumberOfArrays()) {
    this->printErr(
      "Input datasets should have the same number of point data arrays");
    return 0;
  }

  if(input0->GetCellData()->GetNumberOfArrays()
     != input1->GetCellData()->GetNumberOfArrays()) {
    this->printErr(
      "Input datasets should have the same number of cell data arrays");
    return 0;
  }

  output->DeepCopy(input0);
  const auto absdiff
    = [](const double a, const double b) { return std::abs(b - a); };

  for(int i = 0; i < output->GetPointData()->GetNumberOfArrays(); ++i) {
    const auto arrname = output->GetPointData()->GetArrayName(i);
    const auto res = output->GetPointData()->GetArray(i);
    if(res->GetNumberOfComponents() != 1) {
      // skip vector fields
      continue;
    }
    const auto arr0 = input0->GetPointData()->GetArray(arrname);
    const auto arr1 = input1->GetPointData()->GetArray(arrname);
    for(int j = 0; j < res->GetNumberOfTuples(); ++j) {
      res->SetTuple1(j, absdiff(arr0->GetTuple1(j), arr1->GetTuple1(j)));
    }
  }

  for(int i = 0; i < output->GetCellData()->GetNumberOfArrays(); ++i) {
    const auto arrname = output->GetCellData()->GetArrayName(i);
    const auto res = output->GetCellData()->GetArray(i);
    if(res->GetNumberOfComponents() != 1) {
      // skip vector fields
      continue;
    }
    const auto arr0 = input0->GetCellData()->GetArray(arrname);
    const auto arr1 = input1->GetCellData()->GetArray(arrname);
    for(int j = 0; j < res->GetNumberOfTuples(); ++j) {
      res->SetTuple1(j, absdiff(arr0->GetTuple1(j), arr1->GetTuple1(j)));
    }
  }

  return 1;
}
