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

static void setNodes(vtkUnstructuredGrid *const nodes,
                     const std::vector<std::array<float, 3>> &compBaryCoords,
                     const std::vector<int> &compBucketId) {

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
}

static void setArcs(vtkUnstructuredGrid *const arcs,
                    const std::vector<std::set<ttk::SimplexId>> &compArcs,
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

static void extractInputMatrix(ttk::Mapper::Matrix &inputMatrix,
                               std::vector<std::string> &arrayNames,
                               vtkPointData *const pd,
                               const std::string &regexp,
                               const int nThreads,
                               const bool useRegexp) {
  if(pd == nullptr) {
    return;
  }

  if(useRegexp) {
    // select all input columns whose name is matching the regexp
    arrayNames.clear();
    const auto n = pd->GetNumberOfArrays();
    arrayNames.reserve(n);
    for(int i = 0; i < n; ++i) {
      const auto name = pd->GetArrayName(i);

      if(name == nullptr) {
        continue;
      }

      // check regexp
      if(!std::regex_match(name, std::regex(regexp))) {
        continue;
      }

      // sanity check
      const auto arr{pd->GetArray(name)};
      if(arr == nullptr || arr->GetNumberOfComponents() != 1) {
        continue;
      }

      // store matching array names
      arrayNames.emplace_back(name);
    }
  }

  std::vector<vtkDataArray *> arrays(arrayNames.size());
  for(size_t i = 0; i < arrayNames.size(); ++i) {
    arrays[i] = pd->GetArray(arrayNames[i].data());
  }

  inputMatrix.alloc(pd->GetNumberOfTuples(), arrayNames.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < inputMatrix.nRows(); ++i) {
    for(size_t j = 0; j < inputMatrix.nCols(); ++j) {
      const auto array{arrays[j]};
      inputMatrix.get(i, j) = array->GetVariantValue(i).ToDouble();
    }
  }
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

  Matrix inputDistMat{};
  vtkNew<vtkPoints> outputPoints{};
  if(this->ReEmbedMapper) {
    ttk::Timer tm{};
    extractInputMatrix(inputDistMat, this->DistanceMatrix,
                       input->GetPointData(), this->DistanceMatrixRegexp,
                       this->threadNumber_, this->SelectMatrixWithRegexp);
    if(inputDistMat.nRows() != inputDistMat.nCols() || inputDistMat.nCols() == 0
       || inputDistMat.nRows() == 0) {
      this->printErr("Invalid input distance matrix");
      return 0;
    }
    outputPoints->SetNumberOfPoints(input->GetNumberOfPoints());
    outputPoints->GetData()->Fill(0.0);
    auto outSegVTU = vtkUnstructuredGrid::SafeDownCast(outputSegmentation);
    if(outSegVTU != nullptr) {
      outSegVTU->SetPoints(outputPoints);
    }
    this->printMsg("Extracted input distance matrix", 1.0, tm.getElapsedTime(),
                   this->threadNumber_, ttk::debug::LineMode::NEW,
                   ttk::debug::Priority::DETAIL);
  }

  std::vector<std::array<float, 3>> compBaryCoords{};
  std::vector<int> compBucketId{};
  std::vector<std::set<ttk::SimplexId>> compArcs{};

  // calling the base module
  ttkVtkTemplateMacro(
    inputScalarField->GetDataType(), triangulation->getType(),
    this->execute(ttkUtils::GetPointer<int>(bucket),
                  ttkUtils::GetPointer<int>(connComp), compBaryCoords,
                  compBucketId, compArcs,
                  ttkUtils::GetPointer<float>(outputPoints->GetData()),
                  inputDistMat, ttkUtils::GetPointer<VTK_TT>(inputScalarField),
                  *static_cast<TTK_TT *>(triangulation->getData())));

  setNodes(outputNodes, compBaryCoords, compBucketId);
  setArcs(outputArcs, compArcs, outputNodes->GetPoints());

  return 1;
}
