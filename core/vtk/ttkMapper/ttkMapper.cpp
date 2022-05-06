#include <ttkMapper.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <regex>

vtkStandardNewMacro(ttkMapper);

ttkMapper::ttkMapper() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(3);
}

int ttkMapper::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMapper::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

void setNodes(vtkUnstructuredGrid *const nodes,
              const std::vector<std::array<float, 3>> &compBaryCoords,
              const std::vector<std::vector<double>> &highDimBaryCoords,
              const std::vector<std::string> &highDimArrayNames,
              const std::vector<int> &compBucketId,
              const bool computeHighDimBarycenters) {

  vtkNew<vtkIntArray> compId{};
  compId->SetName("ComponentId");
  compId->SetNumberOfComponents(1);
  compId->SetNumberOfTuples(compBaryCoords.size());

  vtkNew<vtkIntArray> bucket{};
  bucket->SetName("BucketId");
  bucket->SetNumberOfComponents(1);
  bucket->SetNumberOfTuples(compBucketId.size());

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(compBaryCoords.size());
  for(size_t i = 0; i < compBaryCoords.size(); ++i) {
    points->SetPoint(i, compBaryCoords[i].data());
    bucket->SetTuple1(i, compBucketId[i]);
    compId->SetTuple1(i, i);
  }

  ttkUtils::CellVertexFromPoints(nodes, points);
  nodes->GetPointData()->AddArray(bucket);
  nodes->GetPointData()->SetScalars(compId);

  if(computeHighDimBarycenters) {
    for(size_t i = 0; i < highDimArrayNames.size(); ++i) {
      vtkNew<vtkDoubleArray> coord{};
      coord->SetName(highDimArrayNames[i].c_str());
      coord->SetNumberOfTuples(highDimBaryCoords.size());
      for(size_t j = 0; j < highDimBaryCoords.size(); ++j) {
        coord->SetTuple1(j, highDimBaryCoords[j][i]);
        nodes->GetPointData()->AddArray(coord);
      }
    }
  }
}

void setArcs(vtkUnstructuredGrid *const arcs,
             const std::vector<std::vector<ttk::SimplexId>> &compArcs,
             vtkPoints *const points) {

  arcs->SetPoints(points);

  std::vector<std::array<vtkIdType, 2>> edges{};
  for(size_t i = 0; i < compArcs.size(); ++i) {
    for(const auto j : compArcs[i]) {
      edges.emplace_back(std::array<vtkIdType, 2>{
        static_cast<vtkIdType>(i), static_cast<vtkIdType>(j)});
    }
  }

  const auto nCells{edges.size()};

  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(nCells + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(2 * nCells);

  vtkNew<vtkIntArray> arcId{};
  arcId->SetName("ArcId");
  arcId->SetNumberOfComponents(1);
  arcId->SetNumberOfTuples(nCells);

  for(size_t i = 0; i < nCells; ++i) {
    offsets->SetTuple1(i, 2 * i);
    connectivity->SetTuple1(2 * i, edges[i][0]);
    connectivity->SetTuple1(2 * i + 1, edges[i][1]);
    arcId->SetTuple1(i, i);
  }
  offsets->SetTuple1(nCells, connectivity->GetNumberOfTuples());

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  arcs->SetCells(VTK_LINE, cells);
  arcs->GetCellData()->SetScalars(arcId);
}

int ttkMapper::RequestData(vtkInformation *ttkNotUsed(request),
                           vtkInformationVector **inputVector,
                           vtkInformationVector *outputVector) {

  if(this->NumberOfBuckets < 1) {
    this->printErr("Invalid requested number of buckets (should be >= 1)");
    return 0;
  }

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputNodes = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputArcs = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputSegmentation = vtkDataSet::GetData(outputVector, 2);

  auto *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    return 0;
  }

  this->preconditionTriangulation(triangulation);

  auto *inputScalarField = this->GetInputArrayToProcess(0, inputVector);
  if(inputScalarField == nullptr) {
    return 0;
  }

  if(inputScalarField->GetNumberOfComponents() != 1) {
    printErr("Invalid scalar field ("
             + std::to_string(inputScalarField->GetNumberOfComponents())
             + " components)");
    return 0;
  }

  const auto pd{input->GetPointData()};
  if(pd != nullptr && this->ComputeHighDimBarycenters
     && this->SelectFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    this->HighDimCoords.clear();
    const auto n = pd->GetNumberOfArrays();
    for(int i = 0; i < n; ++i) {
      const auto &name = pd->GetArrayName(i);
      if(std::regex_match(name, std::regex(RegexpString))) {
        this->HighDimCoords.emplace_back(name);
      }
    }
  }

  std::vector<double> inputHighDimCoords{};
  if(this->ComputeHighDimBarycenters) {
    for(int i = 0; i < input->GetNumberOfPoints(); ++i) {
      for(const auto &arrName : this->HighDimCoords) {
        const auto array{pd->GetArray(arrName.data())};
        inputHighDimCoords.emplace_back(array->GetVariantValue(i).ToDouble());
      }
    }
  }

  vtkNew<vtkIntArray> connComp{};
  connComp->SetName("ComponentId");
  connComp->SetNumberOfComponents(1);
  connComp->SetNumberOfTuples(inputScalarField->GetNumberOfTuples());

  vtkNew<vtkIntArray> bucket{};
  bucket->SetName("BucketId");
  bucket->SetNumberOfComponents(1);
  bucket->SetNumberOfTuples(inputScalarField->GetNumberOfTuples());

  outputSegmentation->ShallowCopy(input);
  outputSegmentation->GetPointData()->AddArray(bucket);
  outputSegmentation->GetPointData()->AddArray(connComp);

  std::vector<std::array<float, 3>> compBaryCoords{};
  std::vector<std::vector<double>> highDimBaryCoords{};
  std::vector<int> compBucketId{};
  std::vector<std::vector<ttk::SimplexId>> compArcs{};

  // calling the base module
  ttkVtkTemplateMacro(
    inputScalarField->GetDataType(), triangulation->getType(),
    this->execute(ttkUtils::GetPointer<int>(bucket),
                  ttkUtils::GetPointer<int>(connComp), compBaryCoords,
                  compBucketId, compArcs, highDimBaryCoords, inputHighDimCoords,
                  ttkUtils::GetPointer<VTK_TT>(inputScalarField),
                  *static_cast<TTK_TT *>(triangulation->getData())));

  setNodes(outputNodes, compBaryCoords, highDimBaryCoords, this->HighDimCoords,
           compBucketId, this->ComputeHighDimBarycenters);
  setArcs(outputArcs, compArcs, outputNodes->GetPoints());

  return 1;
}
