#include <ttkDiscreteMorseFunction.h>

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <vtkUnstructuredGridBase.h>

vtkStandardNewMacro(ttkDiscreteMorseFunction);

ttkDiscreteMorseFunction::ttkDiscreteMorseFunction() {
  this->setDebugMsgPrefix("DiscreteMorseFunction");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkDiscreteMorseFunction::FillInputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkDiscreteMorseFunction::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkDiscreteMorseFunction::RequestData(vtkInformation *ttkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {

  auto *inputDataSet = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto *outputDataSet = vtkUnstructuredGrid::GetData(outputVector, 0);
  if(inputDataSet == nullptr) {
    return 0;
  }

  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, inputVector);
  if(inputScalars == nullptr) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }

  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(inputScalars->GetNumberOfComponents() != 1) {
    this->printErr("Input array needs to be a scalar array.");
    return 0;
  }

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");
  this->printMsg("  Scalar Array: " + std::string(inputScalars->GetName()));

  auto triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(triangulation == nullptr) {
    return 0;
  }
  const auto &explTri
    = *static_cast<ttk::ExplicitTriangulation *>(triangulation->getData());

  this->preconditionTriangulation(triangulation);

  auto inputOffsets = ttkAlgorithm::GetOrderArray(
    inputDataSet, 0, 1, this->ForceInputOffsetScalarField);

  this->setInputScalarField(
    ttkUtils::GetVoidPointer(inputScalars), inputScalars->GetMTime());
  this->setInputOffsets(ttkUtils::GetPointer<SimplexId>(inputOffsets));

  auto status = this->buildGradient(explTri);

  if(status != 0) {
    this->printErr("Could not compute the Discrete Gradient");
    return 0;
  }

  this->dmf = this->computeDiscreteMorseFunction(
    ttkUtils::GetPointer<SimplexId>(inputOffsets), explTri);
  const auto res = this->computeFiltrationOrder(
    ttkUtils::GetPointer<SimplexId>(inputOffsets), explTri);
  this->filtr.resize(res.size());
  for(size_t i = 0; i < res.size(); ++i) {
    this->filtr[res[i].cellId_] = i;
  }

  if(this->dmf.empty()) {
    this->printErr("Could not compute Discrete Morse Function");
    this->dmf.resize(this->filtr.size(), 0);
  }

  if(status != 0) {
    this->printErr("Could not extract the Discrete Morse Function");
    return 0;
  }

  const auto dim = explTri.getDimensionality();
  const auto nVerts = explTri.getNumberOfVertices();
  const auto nEdges = explTri.getNumberOfEdges();
  const auto nTris = dim > 1 ? explTri.getNumberOfTriangles() : 0;
  const auto nTetras = dim > 2 ? explTri.getNumberOfCells() : 0;

  const int nCells = this->dmf.size();
  const int nCellConn = nVerts + 2 * nEdges + 3 * nTris + 4 * nTetras;

  vtkNew<ttkSimplexIdTypeArray> dmfScalars{};
  dmfScalars->SetName(this->DMFArrayName.c_str());
  ttkUtils::SetVoidArray(dmfScalars, this->dmf.data(), this->dmf.size(), 1);

  vtkNew<ttkSimplexIdTypeArray> fltScalars{};
  fltScalars->SetName(this->FiltrArrayName.c_str());
  ttkUtils::SetVoidArray(fltScalars, this->filtr.data(), this->filtr.size(), 1);

  vtkNew<ttkSimplexIdTypeArray> dmfScalarsVerts{};
  dmfScalarsVerts->SetName(this->DMFArrayName.c_str());
  ttkUtils::SetVoidArray(dmfScalarsVerts, this->dmf.data(), nVerts, 1);

  vtkNew<ttkSimplexIdTypeArray> fltScalarsVerts{};
  fltScalarsVerts->SetName(this->FiltrArrayName.c_str());
  ttkUtils::SetVoidArray(fltScalarsVerts, this->filtr.data(), nVerts, 1);

  vtkNew<vtkSignedCharArray> cellDimensions{};
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");
  cellDimensions->SetNumberOfTuples(nCells);

  vtkNew<vtkUnsignedCharArray> cellTypes{};
  cellTypes->SetNumberOfComponents(1);
  cellTypes->SetName("CellTypes");
  cellTypes->SetNumberOfTuples(nCells);

  vtkNew<ttkSimplexIdTypeArray> cellIds{};
  cellIds->SetNumberOfComponents(1);
  cellIds->SetName("CellId");
  cellIds->SetNumberOfTuples(nCells);

  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(nCells + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(nCellConn);
  offsets->SetTuple1(nCells, connectivity->GetNumberOfTuples());

  for(int i = 0; i < nCells; ++i) {
    int o{}, a{};
    std::array<SimplexId, 4> verts{};
    if(i < nVerts) {
      o = i;
      a = o;
      connectivity->SetTuple1(i, i);
      cellDimensions->SetTuple1(i, 0);
      cellTypes->SetTuple1(i, VTK_VERTEX);
    } else if(i < nVerts + nEdges) {
      o = i - nVerts;
      a = 2 * o + nVerts;
      explTri.getEdgeVertex(o, 0, verts[0]);
      explTri.getEdgeVertex(o, 1, verts[1]);
      connectivity->SetTuple1(a + 0, verts[0]);
      connectivity->SetTuple1(a + 1, verts[1]);
      cellDimensions->SetTuple1(i, 1);
      cellTypes->SetTuple1(i, VTK_LINE);
    } else if(i < nVerts + nEdges + nTris) {
      o = i - nVerts - nEdges;
      a = 3 * o + 2 * nEdges + nVerts;
      explTri.getTriangleVertex(o, 0, verts[0]);
      explTri.getTriangleVertex(o, 1, verts[1]);
      explTri.getTriangleVertex(o, 2, verts[2]);
      connectivity->SetTuple1(a + 0, verts[0]);
      connectivity->SetTuple1(a + 1, verts[1]);
      connectivity->SetTuple1(a + 2, verts[2]);
      cellDimensions->SetTuple1(i, 2);
      cellTypes->SetTuple1(i, VTK_TRIANGLE);
    } else {
      o = i - nVerts - nEdges - nTris;
      a = 4 * o + 3 * nTris + 2 * nEdges + nVerts;
      explTri.getCellVertex(o, 0, verts[0]);
      explTri.getCellVertex(o, 1, verts[1]);
      explTri.getCellVertex(o, 2, verts[2]);
      explTri.getCellVertex(o, 3, verts[3]);
      connectivity->SetTuple1(a + 0, verts[0]);
      connectivity->SetTuple1(a + 1, verts[1]);
      connectivity->SetTuple1(a + 2, verts[2]);
      connectivity->SetTuple1(a + 3, verts[3]);
      cellDimensions->SetTuple1(i, 3);
      cellTypes->SetTuple1(i, VTK_TETRA);
    }
    cellIds->SetTuple1(i, o);
    offsets->SetTuple1(i, a);
  }

  // shallow copy input points
  vtkNew<vtkPoints> points{};
  points->ShallowCopy(inputDataSet->GetPoints());
  outputDataSet->SetPoints(points);

  // set point data arrays
  outputDataSet->GetPointData()->AddArray(dmfScalarsVerts);
  outputDataSet->GetPointData()->AddArray(fltScalarsVerts);

  // set cells
  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  outputDataSet->SetCells(cellTypes, cells);

  // set cell data arrays
  outputDataSet->GetCellData()->AddArray(dmfScalars);
  outputDataSet->GetCellData()->AddArray(fltScalars);
  outputDataSet->GetCellData()->AddArray(cellDimensions);
  outputDataSet->GetCellData()->AddArray(cellIds);

  return 1;
}
