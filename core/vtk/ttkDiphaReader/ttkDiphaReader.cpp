#include <ttkDiphaReader.h>
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
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

vtkStandardNewMacro(ttkDiphaReader);

void ttkDiphaReader::PrintSelf(std::ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)")
     << std::endl;
}

ttkDiphaReader::ttkDiphaReader() {
  this->setDebugMsgPrefix("DiphaReader");
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

template <typename T>
void readBin(std::ifstream &stream, T &res) {
  stream.read(reinterpret_cast<char *>(&res), sizeof(res));
}

/**
 * @brief ad-hoc struct to store pairs from a .gudhi file
 */
struct PersistencePair {
  double birth{};
  double death{};
  ttk::SimplexId dim{};
  bool isFinite{true};

  PersistencePair() = default;

  /**
   * @brief comparison operator for sorting pairs
   */
  friend bool operator<(const PersistencePair &lhs,
                        const PersistencePair &rhs) {
    return std::tie(lhs.birth, lhs.death) < std::tie(rhs.birth, rhs.death);
  }

  friend std::ostream &operator<<(std::ostream &os, const PersistencePair &p) {
    os << p.dim << " " << p.birth << " " << p.death << " " << p.isFinite;
    return os;
  }
};

/**
 * @brief diagram is a SORTED vector of pairs (diagonal not included)
 */
using PersistenceDiagram = std::vector<PersistencePair>;

int ttkDiphaReader::readPersistenceDiagram(std::ifstream &stream,
                                           vtkUnstructuredGrid *output) const {

  int64_t nPairs{};
  readBin<int64_t>(stream, nPairs);
  if(nPairs < 0) {
    this->printErr("Negative number of persistence pairs");
    return 0;
  }

  PersistenceDiagram diag(nPairs);

  for(int64_t i = 0; i < nPairs; ++i) {
    int64_t dim{};
    auto &pair = diag[i];
    readBin(stream, dim);
    pair.dim = dim;
    readBin(stream, pair.birth);
    readBin(stream, pair.death);
  }

  // sort pairs
  TTK_PSORT(this->threadNumber_, diag.begin(), diag.end());

  // dataset dimensionality - 1
  const auto rdim = diag.back().dim;

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
  vtkNew<vtkIntArray> critType{};
  critType->SetNumberOfComponents(1);
  critType->SetName("CriticalType");
  critType->SetNumberOfTuples(2 * nPairs);
  vtkNew<vtkDoubleArray> births{};
  births->SetNumberOfComponents(1);
  births->SetName("Birth");
  births->SetNumberOfTuples(nPairs);
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
  vtkNew<vtkUnsignedCharArray> isFinite{};
  isFinite->SetNumberOfComponents(1);
  isFinite->SetName("IsFinite");
  isFinite->SetNumberOfTuples(nPairs);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diag.size(); ++i) {
    const auto &pair = diag[i];
    points->SetPoint(2 * i + 0, pair.birth, pair.birth, 0);
    points->SetPoint(2 * i + 1, pair.birth, pair.death, 0);

    if(pair.dim == -1) { // global extrema pair
      critType->SetTuple1(2 * i + 0, 0);
      critType->SetTuple1(2 * i + 1, 3);
      pairType->SetTuple1(i, 0);
      isFinite->SetTuple1(i, 0);
    } else if(pair.dim < 0) { // essential class
      critType->SetTuple1(2 * i + 0, pair.dim);
      critType->SetTuple1(2 * i + 1, pair.dim + 1);
      pairType->SetTuple1(i, -pair.dim - 1);
      isFinite->SetTuple1(i, 0);
    } else { // regular pairs
      critType->SetTuple1(
        2 * i + 0, (pair.dim == rdim && rdim > 0) ? 2 : pair.dim);
      critType->SetTuple1(2 * i + 1, pair.dim == rdim ? 3 : pair.dim + 1);
      pairType->SetTuple1(i, pair.dim);
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
  const auto maxDeath
    = std::max_element(diag.begin(), diag.end(),
                       [](const PersistencePair &a, const PersistencePair &b) {
                         return a.death < b.death;
                       })
        ->death;
  pairPers->InsertNextTuple1(2 * (maxDeath - minBirth));
  births->InsertNextTuple1(minBirth);

  // copy mesh to output (segfault workaround)
  output->ShallowCopy(mesh);
  // add data arrays
  output->GetPointData()->AddArray(vertId);
  output->GetPointData()->AddArray(critType);
  output->GetPointData()->AddArray(coords);

  output->GetCellData()->AddArray(births);
  output->GetCellData()->AddArray(pairId);
  output->GetCellData()->AddArray(pairType);
  output->GetCellData()->AddArray(pairPers);
  output->GetCellData()->AddArray(isFinite);

  return 1;
}

int ttkDiphaReader::readImageData(std::ifstream &stream,
                                  vtkUnstructuredGrid *output) const {

  // total number of values
  int64_t nVerts{};
  readBin<int64_t>(stream, nVerts);

  // dimension
  int64_t dim{};
  readBin<int64_t>(stream, dim);

  // vertex resolution
  std::array<int64_t, 3> dims{};
  for(int i = 0; i < dim; ++i) {
    readBin<int64_t>(stream, dims[i]);
  }

  // cell resolution
  std::array<int64_t, 3> nCells{
    dims[0] - 1, dims[1] - 1, dim == 3 ? dims[2] - 1 : 1};

  // number of points per cell
  const auto npointscell = dim == 3 ? 8 : 4;

  // grid points
  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(nVerts);
  // grid cells connectivity arrays
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(nCells[0] * nCells[1] * nCells[2] + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(npointscell * nCells[0] * nCells[1]
                                  * nCells[2]);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int64_t k = 0; k < dims[2]; ++k) {
    for(int64_t j = 0; j < dims[1]; ++j) {
      for(int64_t i = 0; i < dims[0]; ++i) {
        points->SetPoint(i + j * dims[0] + k * dims[0] * dims[1], i, j, k);
      }
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int64_t k = 0; k < nCells[2]; ++k) {
    for(int64_t j = 0; j < nCells[1]; ++j) {
      for(int64_t i = 0; i < nCells[0]; ++i) {
        // number of points per grid line
        const auto nptl = dims[0];
        // number of points per grid face
        const auto nptf = nptl * dims[1];
        const auto curr{static_cast<vtkIdType>(i + j * nptl + k * nptf)};
        const auto o = i + j * nCells[0] + k * nCells[0] * nCells[1];
        // build the cell connectivity and offsets array
        if(dim == 3) {
          connectivity->SetTuple1(npointscell * o + 0, curr);
          connectivity->SetTuple1(npointscell * o + 1, curr + 1);
          connectivity->SetTuple1(npointscell * o + 2, curr + nptl);
          connectivity->SetTuple1(npointscell * o + 3, curr + nptl + 1);

          connectivity->SetTuple1(npointscell * o + 4, curr + nptf);
          connectivity->SetTuple1(npointscell * o + 5, curr + nptf + 1);
          connectivity->SetTuple1(npointscell * o + 6, curr + nptf + nptl);
          connectivity->SetTuple1(npointscell * o + 7, curr + nptf + nptl + 1);
        } else if(dim == 2) {
          connectivity->SetTuple1(npointscell * o + 0, curr);
          connectivity->SetTuple1(npointscell * o + 1, curr + 1);
          connectivity->SetTuple1(npointscell * o + 2, curr + nptl);
          connectivity->SetTuple1(npointscell * o + 3, curr + nptl + 1);
        }
        offsets->SetTuple1(o, npointscell * o);
      }
    }
  }
  offsets->SetTuple1(
    nCells[0] * nCells[1] * nCells[2], connectivity->GetNumberOfTuples());

  // scalar field
  vtkNew<vtkDoubleArray> sf{};
  sf->SetName("ScalarField");
  sf->SetNumberOfComponents(1);
  sf->SetNumberOfTuples(nVerts);
  for(int i = 0; i < nVerts; ++i) {
    double val{};
    // ! reading not parallel
    readBin<double>(stream, val);
    sf->SetTuple1(i, val);
  }

  // gather arrays to make the UnstructuredGrid
  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  output->SetPoints(points);
  output->SetCells(dim == 3 ? VTK_VOXEL : VTK_PIXEL, cells);
  output->GetPointData()->AddArray(sf);

  return 1;
}

int ttkDiphaReader::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **ttkNotUsed(inputVector),
                                vtkInformationVector *outputVector) {

  // read data from input file
  std::ifstream stream(FileName, ios::in | ios::binary);

  if(!stream.is_open()) {
    this->printErr("Can't read file: '" + std::string{FileName} + "'");
    return 0;
  }

  int64_t magic{};
  readBin<int64_t>(stream, magic);
  if(magic != 8067171840) {
    this->printErr("Dipha magic number not detected");
    return 0;
  }

  auto output = vtkUnstructuredGrid::GetData(outputVector);

  int64_t type{};
  readBin<int64_t>(stream, type);
  if(type == 1) {
    return readImageData(stream, output);
  } else if(type == 2) {
    return readPersistenceDiagram(stream, output);
  } else if(type == 0) {
    this->printErr("Cannot read Dipha Weighted Boundary Matrices");
    return 0;
  } else if(type == 7 || type == 8) {
    this->printErr("Cannot read Dipha Distance Matrices");
    return 0;
  } else {
    this->printErr("Cannot read unknown Dipha file");
    return 0;
  }
}
