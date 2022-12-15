#include <AbstractTriangulation.h>
#include <RegularGridTriangulation.h>

using namespace std;
using namespace ttk;

AbstractTriangulation::AbstractTriangulation() {

  setDebugMsgPrefix("AbstractTriangulation");

  clear();
}

AbstractTriangulation::~AbstractTriangulation() = default;

int AbstractTriangulation::clear() {

  hasPeriodicBoundaries_ = false;
  hasPreconditionedBoundaryEdges_ = false;
  hasPreconditionedBoundaryTriangles_ = false;
  hasPreconditionedBoundaryVertices_ = false;
  hasPreconditionedCellEdges_ = false;
  hasPreconditionedCellNeighbors_ = false;
  hasPreconditionedCellTriangles_ = false;
  hasPreconditionedEdges_ = false;
  hasPreconditionedEdgeLinks_ = false;
  hasPreconditionedEdgeStars_ = false;
  hasPreconditionedEdgeTriangles_ = false;
  hasPreconditionedTriangles_ = false;
  hasPreconditionedTriangleEdges_ = false;
  hasPreconditionedTriangleLinks_ = false;
  hasPreconditionedTriangleStars_ = false;
  hasPreconditionedVertexEdges_ = false;
  hasPreconditionedVertexLinks_ = false;
  hasPreconditionedVertexNeighbors_ = false;
  hasPreconditionedVertexStars_ = false;
  hasPreconditionedVertexTriangles_ = false;
  hasPreconditionedManifold_ = false;

  boundaryEdges_.clear();
  boundaryTriangles_.clear();
  boundaryVertices_.clear();

  tetraEdgeList_.clear();
  cellNeighborList_.clear();
  tetraTriangleList_.clear();

  edgeLinkList_.clear();
  edgeList_.clear();
  edgeStarList_.clear();
  edgeTriangleList_.clear();

  triangleList_.clear();
  triangleEdgeList_.clear();
  triangleLinkList_.clear();
  triangleStarList_.clear();

  vertexEdgeList_.clear();
  vertexLinkList_.clear();
  vertexNeighborList_.clear();
  vertexStarList_.clear();
  vertexTriangleList_.clear();

  return 0;
}

template <class itemType>
size_t AbstractTriangulation::tableTableFootprint(
  const vector<vector<itemType>> &table,
  const string &tableName,
  ostream &stream) const {

  size_t localByteNumber = 0;
  stringstream msg;

  for(size_t i = 0; i < table.size(); i++) {
    localByteNumber += table[i].size() * sizeof(itemType);
  }

  if((localByteNumber) && (tableName.length()) && (msg)) {
    msg << tableName << ": " << localByteNumber << " bytes";
    printMsg(msg.str(), debug::Priority::INFO, debug::LineMode::NEW, stream);
  }

  return localByteNumber;
}

size_t AbstractTriangulation::footprint(size_t size) const {

  size += sizeof(*this);
  stringstream msg;

  size += tableFootprint<bool>(boundaryEdges_, "boundaryEdges_");

  size += tableFootprint<bool>(boundaryTriangles_, "boundaryTriangles_");

  size += tableFootprint<bool>(boundaryVertices_, "boundaryVertices_");

  size += tableFootprint(tetraEdgeList_, "tetraEdgeList_");

  size
    += tableTableFootprint<SimplexId>(cellNeighborList_, "cellNeighborList_");

  size += tableFootprint(tetraTriangleList_, "tetraTriangleList_");

  size += tableTableFootprint<SimplexId>(edgeLinkList_, "edgeLinkList_");

  size += tableFootprint(edgeList_, "edgeList_");

  size += tableTableFootprint<SimplexId>(edgeStarList_, "edgeStarList_");

  size
    += tableTableFootprint<SimplexId>(edgeTriangleList_, "edgeTriangleList_");

  size += tableFootprint(triangleList_, "triangleList_");

  size += tableFootprint(triangleEdgeList_, "triangleEdgeList_");

  size
    += tableTableFootprint<SimplexId>(triangleLinkList_, "triangleLinkList_");

  size
    += tableTableFootprint<SimplexId>(triangleStarList_, "triangleStarList_");

  size += tableTableFootprint<SimplexId>(vertexEdgeList_, "vertexEdgeList_");

  size += tableTableFootprint<SimplexId>(vertexLinkList_, "vertexLinkList_");

  size += tableTableFootprint<SimplexId>(
    vertexNeighborList_, "vertexNeighborList_");

  size += tableTableFootprint<SimplexId>(vertexStarList_, "vertexStarList_");

  size += tableTableFootprint<SimplexId>(
    vertexTriangleList_, "vertexTriangleList_");

  size += tableTableFootprint(cellEdgeVector_, "cellEdgeVector_");
  size += tableTableFootprint(cellTriangleVector_, "cellTriangleVector_");
  size += tableTableFootprint(triangleEdgeVector_, "triangleEdgeVector_");

  msg << "Total footprint: " << (size / 1024) / 1024 << " MB.";
  printMsg(msg.str());

  return size;
}

SimplexId ttk::RegularGridTriangulation::findEdgeFromVertices(
  const SimplexId v0, const SimplexId v1) const {
  // loop over v0 edges to find the one between v0 and v1
  const auto nEdges = this->getVertexEdgeNumberInternal(v0);
  for(SimplexId i = 0; i < nEdges; ++i) {
    SimplexId e{};
    std::array<SimplexId, 2> eVerts{};
    this->getVertexEdgeInternal(v0, i, e);
    this->getEdgeVertexInternal(e, 0, eVerts[0]);
    this->getEdgeVertexInternal(e, 1, eVerts[1]);
    if((v0 == eVerts[0] && v1 == eVerts[1])
       || (v0 == eVerts[1] && v1 == eVerts[0])) {
      return e;
    }
  }

  return -1;
}

SimplexId ttk::RegularGridTriangulation::findTriangleFromVertices(
  std::array<SimplexId, 3> &verts) const {

  std::sort(verts.begin(), verts.end());

  // loop over verts[0] triangles to find the one shared by all 3
  const auto nTriangles = this->getVertexTriangleNumberInternal(verts[0]);
  for(SimplexId i = 0; i < nTriangles; ++i) {
    SimplexId t{};
    std::array<SimplexId, 3> tVerts{};
    this->getVertexTriangleInternal(verts[0], i, t);
    this->getTriangleVertexInternal(t, 0, tVerts[0]);
    this->getTriangleVertexInternal(t, 1, tVerts[1]);
    this->getTriangleVertexInternal(t, 2, tVerts[2]);
    std::sort(tVerts.begin(), tVerts.end());
    if(tVerts == verts) {
      return t;
    }
  }

  return -1;
}
