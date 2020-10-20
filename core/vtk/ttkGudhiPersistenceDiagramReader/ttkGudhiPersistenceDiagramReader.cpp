#include <ttkGudhiPersistenceDiagramReader.h>
#include <ttkMacros.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkLongLongArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

vtkStandardNewMacro(ttkGudhiPersistenceDiagramReader);

void ttkGudhiPersistenceDiagramReader::PrintSelf(std::ostream &os,
                                                 vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)")
     << std::endl;
}

ttkGudhiPersistenceDiagramReader::ttkGudhiPersistenceDiagramReader() {
  this->setDebugMsgPrefix("GudhiPersistenceDiagramReader");
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

// ad-hoc struct to store pairs from a .gudhi file
struct PersistencePair {
  double birth{};
  double death{};
  ttk::SimplexId dim{};
  bool isFinite{true};

  PersistencePair(const double bi,
                  const double de,
                  const ttk::SimplexId di,
                  const bool fi)
    : birth(bi), death(de), dim(di), isFinite(fi) {
  }

  // comparison operator for sorting pairs
  friend bool operator<(const PersistencePair &lhs,
                        const PersistencePair &rhs) {
    return std::tie(lhs.birth, lhs.death) < std::tie(rhs.birth, rhs.death);
  }

  friend std::ostream &operator<<(std::ostream &os, const PersistencePair &p) {
    os << p.dim << " " << p.birth << " " << p.death << " " << p.isFinite;
    return os;
  }
};

// diagram is a SORTED vector of pairs (diagonal not included)
using PersistenceDiagram = std::vector<PersistencePair>;

int ttkGudhiPersistenceDiagramReader::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **ttkNotUsed(inputVector),
  vtkInformationVector *outputVector) {

  // read data from input file
  std::ifstream stream(FileName, ios::in);

  if(!stream.is_open()) {
    this->printErr("Can't read file: '" + std::string{FileName} + "'");
    return -1;
  }

  PersistenceDiagram diag{};

  std::string line{};
  while(!stream.eof()) {
    std::getline(stream, line);

    // ignore empty and commented lines (beginning with a "#")
    if(line.empty() || line[0] == '#') {
      continue;
    }

    std::array<double, 4> values{};
    std::stringstream ss(line);
    int nVals = std::sscanf(line.c_str(), "%lf %lf %lf %lf", &values[0],
                            &values[1], &values[2], &values[3]);
    if(nVals > 1) {
      const auto isFinite = std::isfinite(values[nVals - 1]);
      diag.emplace_back(values[nVals - 2], isFinite ? values[nVals - 1] : -1,
                        nVals > 2 ? values[nVals - 3] : -1, isFinite);
    }
  }

  if(diag.empty()) {
    this->printErr("Diagram has no pairs");
    return 0;
  }

  const auto nPairs = diag.size();

  // convert data to UnstructuredGrid
  vtkNew<vtkUnstructuredGrid> mesh{};
  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(2 * nPairs);
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(nPairs + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(2 * nPairs);

  vtkNew<vtkIntArray> pairType{};
  pairType->SetNumberOfComponents(1);
  pairType->SetName("PairType");
  pairType->SetNumberOfTuples(nPairs);
  vtkNew<ttkSimplexIdTypeArray> pairId{};
  pairId->SetNumberOfComponents(1);
  pairId->SetName("PairIdentifier");
  pairId->SetNumberOfTuples(nPairs);
  vtkNew<vtkDoubleArray> pairPers{};
  pairPers->SetNumberOfComponents(1);
  pairPers->SetName("Persistence");
  pairPers->SetNumberOfTuples(nPairs);
  vtkNew<vtkUnsignedCharArray> isFinite{};
  isFinite->SetNumberOfComponents(1);
  isFinite->SetName("IsFinite");
  isFinite->SetNumberOfTuples(nPairs);
  vtkNew<vtkDoubleArray> births{};
  births->SetNumberOfComponents(1);
  births->SetName("Birth");
  births->SetNumberOfTuples(nPairs);

  vtkNew<vtkIntArray> critType{};
  critType->SetNumberOfComponents(1);
  critType->SetName("CriticalType");
  critType->SetNumberOfTuples(2 * nPairs);
  // misc dummy point data
  vtkNew<ttkSimplexIdTypeArray> vertId{};
  vertId->SetNumberOfComponents(1);
  vertId->SetName(ttk::VertexScalarFieldName);
  vertId->SetNumberOfTuples(2 * nPairs);
  vertId->Fill(0);
  vtkNew<vtkFloatArray> coords{};
  coords->SetNumberOfComponents(3);
  coords->SetName("Coordinates");
  coords->SetNumberOfTuples(2 * nPairs);
  coords->Fill(0.0F);

  const auto maxDeath
    = std::max_element(diag.begin(), diag.end(),
                       [](const PersistencePair &a, const PersistencePair &b) {
                         return a.death < b.death;
                       })
        ->death;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diag.size(); ++i) {
    auto &pair = diag[i];
    if(pair.dim == -1 || pair.death == -1) {
      pair.isFinite = false;
      pair.death = maxDeath;
    }
  }

  // sort pairs
  TTK_PSORT(this->threadNumber_, diag.begin(), diag.end());

  // dataset dimensionality - 1
  const auto rdim = diag.back().dim;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nPairs; ++i) {
    const auto &pair = diag[i];
    points->SetPoint(2 * i + 0, pair.birth, pair.birth, 0);
    points->SetPoint(2 * i + 1, pair.birth, pair.death, 0);

    critType->SetTuple1(
      2 * i + 0, (pair.dim == rdim && rdim > 0) ? 2 : pair.dim);
    critType->SetTuple1(2 * i + 1, pair.dim == rdim ? 3 : pair.dim + 1);
    pairType->SetTuple1(i, pair.dim);

    if(!pair.isFinite) { // global extrema pair
      isFinite->SetTuple1(i, 0);
    } else { // regular pairs
      isFinite->SetTuple1(i, 1);
    }

    connectivity->SetTuple1(2 * i, 2 * i);
    connectivity->SetTuple1(2 * i + 1, 2 * i + 1);
    offsets->SetTuple1(i, 2 * i);
    pairId->SetTuple1(i, i);
    pairPers->SetTuple1(i, pair.death - pair.birth);
    births->SetTuple1(i, pair.birth);
  }

  offsets->SetTuple1(nPairs, connectivity->GetNumberOfTuples());
  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  mesh->SetCells(VTK_LINE, cells);
  mesh->SetPoints(points);

  // diagonal
  const auto minBirth = diag[0].birth;
  std::array<vtkIdType, 2> diagIds{
    // id of global minimum in vtkPoints
    0,
    // id of local minimum with highest birth in vtkPoints
    static_cast<vtkIdType>(2 * (nPairs - 1)),
  };
  mesh->InsertNextCell(VTK_LINE, 2, diagIds.data());
  pairId->InsertNextTuple1(-1); // diagonal id = -1
  pairType->InsertNextTuple1(-1);
  isFinite->InsertNextTuple1(0);
  pairPers->InsertNextTuple1(2 * (maxDeath - minBirth));
  births->InsertNextTuple1(minBirth);

  // copy mesh to output (segfault workaround)
  auto output = vtkUnstructuredGrid::GetData(outputVector);
  output->ShallowCopy(mesh);
  // add data arrays
  output->GetPointData()->AddArray(vertId);
  output->GetPointData()->AddArray(critType);
  output->GetPointData()->AddArray(coords);

  output->GetCellData()->AddArray(pairId);
  output->GetCellData()->AddArray(pairType);
  output->GetCellData()->AddArray(pairPers);
  output->GetCellData()->AddArray(births);
  output->GetCellData()->AddArray(isFinite);

  return 1;
}
