#include <AbstractTriangulation.h>

using namespace ttk;

AbstractTriangulation::AbstractTriangulation() {
  setDebugMsgPrefix("AbstractTriangulation");
  clear();
}

AbstractTriangulation::~AbstractTriangulation() = default;

void AbstractTriangulation::clear() {
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
}

size_t AbstractTriangulation::footprint(size_t size) const {

  size += sizeof(*this);

  const auto printArrayFootprint
    = [this](const FlatJaggedArray &array, const std::string &name) {
        if(!array.empty() && !name.empty()) {
          this->printMsg(name + std::string{": "}
                         + std::to_string(array.footprint()) + " bytes");
        }
        return array.footprint();
      };

  size += printArrayFootprint(vertexNeighborList_, "vertexNeighborList_");
  size += printArrayFootprint(cellNeighborList_, "cellNeighborList_");
  size += printArrayFootprint(vertexEdgeList_, "vertexEdgeList_");
  size += printArrayFootprint(vertexTriangleList_, "vertexTriangleList_");
  size += printArrayFootprint(edgeTriangleList_, "edgeTriangleList_");
  size += printArrayFootprint(vertexStarList_, "vertexStarList_");
  size += printArrayFootprint(edgeStarList_, "edgeStarList_");
  size += printArrayFootprint(triangleStarList_, "triangleStarList_");
  size += printArrayFootprint(vertexLinkList_, "vertexLinkList_");
  size += printArrayFootprint(edgeLinkList_, "edgeLinkList_");
  size += printArrayFootprint(triangleLinkList_, "triangleLinkList_");

  size += tableFootprint(edgeList_, "edgeList_");
  size += tableFootprint(triangleList_, "triangleList_");
  size += tableFootprint(triangleEdgeList_, "triangleEdgeList_");
  size += tableFootprint(tetraEdgeList_, "tetraEdgeList_");
  size += tableFootprint(tetraTriangleList_, "tetraTriangleList_");

  size += tableFootprint(boundaryVertices_, "boundaryVertices_");
  size += tableFootprint(boundaryEdges_, "boundaryEdges_");
  size += tableFootprint(boundaryTriangles_, "boundaryTriangles_");

  this->printMsg("Total footprint: " + std::to_string((size / 1024) / 1024)
                 + " MiB.");

  return size;
}
