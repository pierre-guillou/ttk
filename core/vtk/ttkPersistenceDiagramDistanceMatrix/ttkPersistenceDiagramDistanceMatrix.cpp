#include <ttkPersistenceDiagramDistanceMatrix.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkTable.h>

vtkStandardNewMacro(ttkPersistenceDiagramDistanceMatrix);

ttkPersistenceDiagramDistanceMatrix::ttkPersistenceDiagramDistanceMatrix() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkPersistenceDiagramDistanceMatrix::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagramDistanceMatrix::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramDistanceMatrix::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Memory m;

  // Get input data
  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  auto nBlocks = inputVector[0]->GetNumberOfInformationObjects();
  std::vector<vtkMultiBlockDataSet *> blocks(nBlocks);

  if(nBlocks > 2) {
    this->printWrn("Only dealing with the first two MultiBlockDataSets");
    nBlocks = 2;
  }

  // number of diagrams per input block
  std::array<size_t, 2> nInputs{0, 0};

  for(int i = 0; i < nBlocks; ++i) {
    blocks[i] = vtkMultiBlockDataSet::GetData(inputVector[0], i);
    if(blocks[i] != nullptr) {
      nInputs[i] = blocks[i]->GetNumberOfBlocks();
      for(size_t j = 0; j < nInputs[i]; ++j) {
        inputDiagrams.emplace_back(
          vtkUnstructuredGrid::SafeDownCast(blocks[i]->GetBlock(j)));
      }
    }
  }

  // total number of diagrams
  const int nDiags = inputDiagrams.size();

  // Sanity check
  for(const auto vtu : inputDiagrams) {
    if(vtu == nullptr) {
      this->printErr("Input diagrams are not all vtkUnstructuredGrid");
      return 0;
    }
  }

  // Set output
  auto diagramsDistTable = vtkTable::GetData(outputVector);

  std::vector<ttk::Diagram> intermediateDiagrams(nDiags);

  for(int i = 0; i < nDiags; i++) {
    const auto ret = VTUToDiagram(intermediateDiagrams[i], inputDiagrams[i]);
    if(ret < 0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
  }

  const auto diagramsDistMat = this->execute(intermediateDiagrams, nInputs);

  // zero-padd column name to keep Row Data columns ordered
  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  const auto nTuples = nInputs[1] == 0 ? nInputs[0] : nInputs[1];

  // copy diagrams distance matrix to output
  for(size_t i = 0; i < diagramsDistMat.size(); ++i) {
    std::string name{"Diagram"};
    zeroPad(name, diagramsDistMat.size(), i);

    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(nTuples);
    col->SetName(name.c_str());
    for(size_t j = 0; j < diagramsDistMat[i].size(); ++j) {
      col->SetTuple1(j, diagramsDistMat[i][j]);
    }
    diagramsDistTable->AddColumn(col);
  }

  // aggregate input field data
  vtkNew<vtkFieldData> fd{};
  fd->CopyStructure(inputDiagrams[0]->GetFieldData());
  fd->SetNumberOfTuples(nTuples);
  for(size_t i = 0; i < nTuples; ++i) {
    fd->SetTuple(i, 0, inputDiagrams[i]->GetFieldData());
  }

  // copy input field data to output row data
  for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
    diagramsDistTable->AddColumn(fd->GetAbstractArray(i));
  }

  return 1;
}

int ttkPersistenceDiagramDistanceMatrix::VTUToDiagram(
  ttk::Diagram &diagram, vtkUnstructuredGrid *vtu) const {

  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();
  const auto points = vtu->GetPoints();

  if(pd == nullptr) {
    this->printErr("VTU diagram with NULL Point Data");
    return -1;
  }
  if(cd == nullptr) {
    this->printErr("VTU diagram with NULL Cell Data");
    return -2;
  }
  if(points == nullptr) {
    this->printErr("VTU with no points");
    return -3;
  }

  // cell data
  const auto pairId = vtkIntArray::SafeDownCast(cd->GetArray("PairIdentifier"));
  const auto pairType = vtkIntArray::SafeDownCast(cd->GetArray("PairType"));
  const auto pairPers
    = vtkDoubleArray::SafeDownCast(cd->GetArray("Persistence"));

  // point data
  const auto vertexId
    = vtkIntArray::SafeDownCast(pd->GetArray(ttk::VertexScalarFieldName));
  const auto critType = vtkIntArray::SafeDownCast(pd->GetArray("CriticalType"));
  const auto birthScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Birth"));
  const auto deathScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Death"));
  const auto coords = vtkFloatArray::SafeDownCast(pd->GetArray("Coordinates"));

  const bool embed = coords == nullptr;

  int nPairs = pairId->GetNumberOfTuples();

  // compact pairIds in [0, nPairs - 1] (diagonal excepted)
  for(int i = 0; i < nPairs; i++) {
    if(pairId->GetTuple1(i) != -1) {
      pairId->SetTuple1(i, i);
    }
  }

  // skip diagram diagonal if present (assuming it's the last pair in the
  // diagram)
  if(pairId->GetTuple1(nPairs - 1) == -1) {
    nPairs -= 1;
  }

  if(nPairs < 1 || vertexId == nullptr || pairId == nullptr
     || critType == nullptr || pairPers == nullptr || pairType == nullptr
     || birthScalars == nullptr || deathScalars == nullptr) {
    this->printErr("Either no pairs in diagram or some array is NULL");
    return -4;
  }

  diagram.resize(nPairs + 1);

  // count the number of pairs whose index is >= nPairs
  int nbNonCompact = 0;

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < nPairs; ++i) {

    const auto i0 = 2 * i + 0;
    const auto i1 = 2 * i + 1;

    const auto v0 = vertexId->GetValue(i0);
    const auto v1 = vertexId->GetValue(i1);
    const auto ct0 = static_cast<ttk::CriticalType>(critType->GetValue(i0));
    const auto ct1 = static_cast<ttk::CriticalType>(critType->GetValue(i1));

    const auto pId = pairId->GetValue(i);
    const auto pType = pairType->GetValue(i);
    const auto pers = pairPers->GetValue(i);

    const auto birth = birthScalars->GetValue(i0);
    const auto death = deathScalars->GetValue(i1);

    std::array<double, 3> coordsBirth{}, coordsDeath{};

    if(embed) {
      points->GetPoint(i0, coordsBirth.data());
      points->GetPoint(i1, coordsDeath.data());
    } else {
      coords->GetTuple(i0, coordsBirth.data());
      coords->GetTuple(i1, coordsDeath.data());
    }

    if(pId != -1 && pId < nPairs) {

      if(pId == 0) {
        // duplicate the global min-max pair into two: one min-saddle pair and
        // one saddle-max pair
        diagram[0] = std::make_tuple(
          v0, ttk::CriticalType::Local_minimum, v1, ttk::CriticalType::Saddle1,
          pers, pType, birth, coordsBirth[0], coordsBirth[1], coordsBirth[2],
          death, coordsDeath[0], coordsDeath[1], coordsDeath[2]);
        // store the saddle max pair at the vector end
        diagram[nPairs] = std::make_tuple(
          v0, ttk::CriticalType::Saddle1, v1, ttk::CriticalType::Local_maximum,
          pers, pType, birth, coordsBirth[0], coordsBirth[1], coordsBirth[2],
          death, coordsDeath[0], coordsDeath[1], coordsDeath[2]);

      } else {
        // all other pairs
        diagram[pId] = std::make_tuple(v0, ct0, v1, ct1, pers, pType, birth,
                                       coordsBirth[0], coordsBirth[1],
                                       coordsBirth[2], death, coordsDeath[0],
                                       coordsDeath[1], coordsDeath[2]);
      }
    }

    if(pId >= nPairs) {
      nbNonCompact++;
    }
  }

  if(nbNonCompact > 0) {
    this->printWrn("Missed " + std::to_string(nbNonCompact)
                   + " pairs due to non-compactness.");
  }

  return 0;
}
