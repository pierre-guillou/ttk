/// \ingroup base
/// \class ttk::Mapper
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date May 2022.
///
/// \brief TTK processing package for mapper.
///
/// This class generates a mapper from a data-set with a given number
/// of buckets. Data-set vertices and edges are placed into buckets
/// from their input scalar field value then connected components are
/// detected per bucket.
///
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkMapper.cpp %for a usage example.

#pragma once

#include <Triangulation.h>
#include <UnionFind.h>

#include <array>

namespace ttk {

  class Mapper : virtual public Debug {
  public:
    Mapper();

    inline void
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation != nullptr) {
        triangulation->preconditionEdges();
        triangulation->preconditionVertexEdges();
      }
    }

    /**
     * @brief Compute mapper
     *
     * @param[out] outputBucket Bucket id for each vertex
     * @param[out] outputConnComps Component id for each vertex
     * @param[out] compBaryCoords 3D barycenters coordinates
     * @param[out] compArcs Arcs between connected components
     * @param[out] highDimBaryCoords (Optional) High-dimensional barycenters
     * coordinates
     * @param[in] inputHighDimCoords (Optional) High-dimensional input points
     * coordinates
     * @param[in] inputSf Input scalar field (on vertices)
     * @param[in] triangulation Triangulation
     *
     * @return 0 in case of success.
     */
    template <typename dataType, typename triangulationType>
    int execute(int *const outputBucket,
                int *const outputConnComp,
                std::vector<std::array<float, 3>> &compBaryCoords,
                std::vector<int> &connCompBucket,
                std::vector<std::vector<SimplexId>> &compArcs,
                std::vector<std::vector<double>> &highDimBaryCoords,
                const std::vector<double> &inputHighDimCoords,
                const dataType *const inputSf,
                const triangulationType &triangulation) const;

  private:
    /**
     * @brief Map vertices & edges to buckets
     *
     * @param[out] vertsBucket Bucket id for each vertex
     * @param[out] bucketEdges List of edges in bucket per bucket
     * @param[in] inputSf Input scalar field (on vertices)
     * @param[in] triangulation Triangulation
     *
     * @return 0 in case of success.
     */
    template <typename dataType, typename triangulationType>
    int findBuckets(int *const vertsBucket,
                    std::vector<std::vector<SimplexId>> &bucketEdges,
                    const dataType *const inputSf,
                    const triangulationType &triangulation) const;

    /**
     * @brief Detect connected components for a given bucket.
     *
     * @param[in] bucketId Id of current bucket
     * @param[out] connComps Component id for each vertex
     * @param[in] vertsBucket Bucket id for each vertex
     * @param[out] connCompEdges List of edges per connected component
     * @param[in] bucketEdges List of edges in the current bucket (reused
     * storage)
     * @param[in] edgeGidToLid Global to compact edge mapping (reused storage)
     * @param[in] uf Local edge Union-Find vector (reused storage)
     * @param[out] connCompIds Local edge connected component id (reused
     * storage)
     * @param[in] triangulation Triangulation
     *
     * @return 0 in case of success.
     */
    template <typename triangulationType>
    int processBucket(const int bucketId,
                      int *const connComps,
                      const int *const vertsBucket,
                      std::vector<std::vector<SimplexId>> &connCompEdges,
                      std::vector<SimplexId> &bucketEdges,
                      std::vector<SimplexId> &edgeGidToLid,
                      std::vector<ttk::UnionFind> &uf,
                      std::vector<int> &connCompIds,
                      const triangulationType &triangulation) const;

    /**
     * @brief Compute connected component barycenter.
     *
     * @param[out] baryCoords 3D coordinates for current component
     * @param[out] highDimBaryCoords (Optional) High-dimensional barycenter
     * coordinates
     * @param[in] inputHighDimCoords (Optional) High-dimensional input points
     * coordinates
     * @param[in] connCompEdges List of edges of the current component
     * @param[in] triangulation Triangulation
     *
     * @return 0 in case of success.
     */
    template <typename triangulationType>
    int computeCompBarycenter(std::array<float, 3> &baryCoords,
                              std::vector<double> &highDimBaryCoords,
                              const std::vector<double> &inputHighDimCoords,
                              const std::vector<SimplexId> &connCompEdges,
                              const triangulationType &triangulation) const;

  protected:
    int NumberOfBuckets{10};
    bool ComputeHighDimBarycenters{false};
  };

} // namespace ttk

// template functions
template <typename dataType, typename triangulationType>
int ttk::Mapper::execute(int *const outputBucket,
                         int *const outputConnComp,
                         std::vector<std::array<float, 3>> &compBaryCoords,
                         std::vector<int> &connCompBucket,
                         std::vector<std::vector<SimplexId>> &compArcs,
                         std::vector<std::vector<double>> &highDimBaryCoords,
                         const std::vector<double> &inputHighDimCoords,
                         const dataType *const inputSf,
                         const triangulationType &triangulation) const {

  Timer tm{}, tmsec{};

  std::vector<std::vector<SimplexId>> bucketEdges(this->NumberOfBuckets);

  this->findBuckets(outputBucket, bucketEdges, inputSf, triangulation);

  // allocate memory
  compBaryCoords.reserve(2 * this->NumberOfBuckets);
  compArcs.reserve(2 * this->NumberOfBuckets);

  std::vector<SimplexId> edgeGidToLid(triangulation.getNumberOfEdges(), -1);
  std::vector<ttk::UnionFind> uf{};
  std::vector<int> connCompIds{};
  std::vector<std::vector<std::vector<SimplexId>>> connCompsEdgesPerBucket(
    this->NumberOfBuckets);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) \
  firstprivate(edgeGidToLid, uf, connCompIds)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < this->NumberOfBuckets; ++i) {
    this->processBucket(i, outputConnComp, outputBucket,
                        connCompsEdgesPerBucket[i], bucketEdges[i],
                        edgeGidToLid, uf, connCompIds, triangulation);
  }

  this->printMsg("Detected connected components", 1.0, tmsec.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
  tmsec.reStart();

  size_t nConnComps{};
  // (bucket id, local component id) -> global component id
  std::map<std::pair<int, size_t>, int> connCompLidToGid{};
  for(int i = 0; i < this->NumberOfBuckets; ++i) {
    for(size_t j = 0; j < connCompsEdgesPerBucket[i].size(); ++j) {
      connCompLidToGid[{i, j}] = nConnComps++;
    }
  }

  // offset component id
  const auto nVerts{triangulation.getNumberOfVertices()};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    outputConnComp[i] = connCompLidToGid[{outputBucket[i], outputConnComp[i]}];
  }

  this->printMsg("Generated segmentation", 1.0, tmsec.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
  tmsec.reStart();

  std::vector<std::vector<SimplexId>> connCompEdges(nConnComps);
  // aggregate component edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < this->NumberOfBuckets; ++i) {
    for(size_t j = 0; j < connCompsEdgesPerBucket[i].size(); ++j) {
      connCompEdges[connCompLidToGid[{i, j}]]
        = std::move(connCompsEdgesPerBucket[i][j]);
    }
  }

  // find barycenter coordinates
  compBaryCoords.resize(nConnComps);
  highDimBaryCoords.resize(nConnComps);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nConnComps; ++i) {
    this->computeCompBarycenter(compBaryCoords[i], highDimBaryCoords[i],
                                inputHighDimCoords, connCompEdges[i],
                                triangulation);
  }

  this->printMsg("Found barycenter coordinates", 1.0, tmsec.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
  tmsec.reStart();

  // transposition of connCompEdges
  std::vector<std::vector<int>> edgeConnComps(triangulation.getNumberOfEdges());
  for(size_t i = 0; i < connCompEdges.size(); ++i) {
    for(const auto e : connCompEdges[i]) {
      edgeConnComps[e].emplace_back(i);
    }
  }

  // global component id -> bucket id
  connCompBucket.resize(nConnComps);
  for(const auto &p : connCompLidToGid) {
    connCompBucket[p.second] = p.first.first;
  }

  // detect arcs between connected components
  compArcs.resize(nConnComps);
  for(size_t i = 0; i < edgeConnComps.size(); ++i) {
    // an arc edge should at least cover two connected components
    if(edgeConnComps[i].size() < 2) {
      continue;
    }
    for(size_t j = 1; j < edgeConnComps[i].size(); ++j) {
      const auto currCompId{edgeConnComps[i][j]};
      for(size_t k = 0; k < j; ++k) {
        const auto prevCompId{edgeConnComps[i][k]};
        if(connCompBucket[prevCompId] == connCompBucket[currCompId] - 1) {
          compArcs[currCompId].emplace_back(prevCompId);
        }
      }
    }
  }

  // remove duplicate arcs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < compArcs.size(); ++i) {
    std::sort(compArcs[i].begin(), compArcs[i].end());
    const auto last = std::unique(compArcs[i].begin(), compArcs[i].end());
    compArcs[i].erase(last, compArcs[i].end());
  }

  this->printMsg("Detected arcs between connected components", 1.0,
                 tmsec.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);
  tmsec.reStart();

  this->printMsg(
    "Found " + std::to_string(nConnComps) + " connected components", 1,
    tm.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::Mapper::findBuckets(int *const vertsBucket,
                             std::vector<std::vector<SimplexId>> &bucketEdges,
                             const dataType *const inputSf,
                             const triangulationType &triangulation) const {

  Timer tm{};

  const auto nVerts{triangulation.getNumberOfVertices()};
  const auto sfRange = std::minmax_element(inputSf, inputSf + nVerts);

  std::vector<double> buckets(this->NumberOfBuckets + 1);
  buckets.front() = *sfRange.first;
  buckets.back() = *sfRange.second;
  const double bucketLength{(buckets.back() - buckets.front())
                            / this->NumberOfBuckets};
  for(int i = 1; i < this->NumberOfBuckets; ++i) {
    buckets[i] = buckets[i - 1] + bucketLength;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    for(int j = 0; j < this->NumberOfBuckets; ++j) {
      if(inputSf[i] < buckets[j + 1]) {
        vertsBucket[i] = j;
        break;
      }
      vertsBucket[i] = this->NumberOfBuckets - 1;
    }
  }

  const auto nEdges{triangulation.getNumberOfEdges()};

  // find edges covering bucket
  for(SimplexId e = 0; e < nEdges; ++e) {
    SimplexId v0{}, v1{};
    triangulation.getEdgeVertex(e, 0, v0);
    triangulation.getEdgeVertex(e, 1, v1);
    auto b0 = vertsBucket[v0];
    auto b1 = vertsBucket[v1];
    if(b0 > b1) {
      std::swap(b0, b1);
    }
    for(int i = b0; i < b1 + 1; ++i) {
      bucketEdges[i].emplace_back(e);
    }
  }

  this->printMsg("Placed vertices & edges in buckets", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
}

template <typename triangulationType>
int ttk::Mapper::processBucket(
  const int bucketId,
  int *const connComps,
  const int *const vertsBucket,
  std::vector<std::vector<SimplexId>> &connCompEdges,
  std::vector<SimplexId> &bucketEdges,
  std::vector<SimplexId> &edgeGidToLid,
  std::vector<ttk::UnionFind> &uf,
  std::vector<int> &connCompIds,
  const triangulationType &triangulation) const {

  // init uf, connCompIds vectors
  uf.resize(bucketEdges.size());
  connCompIds.resize(bucketEdges.size());
  for(size_t i = 0; i < bucketEdges.size(); ++i) {
    edgeGidToLid[bucketEdges[i]] = i;
    uf[i].setParent(&uf[i]);
    uf[i].setRank(i);
    connCompIds[i] = -1;
  }

  // find connected components in edges
  for(size_t i = 0; i < bucketEdges.size(); ++i) {
    const auto e{bucketEdges[i]};
    for(SimplexId j = 0; j < 2; ++j) {
      SimplexId v{};
      triangulation.getEdgeVertex(e, j, v);
      for(SimplexId k = 0; k < triangulation.getVertexEdgeNumber(v); ++k) {
        SimplexId oe{};
        triangulation.getVertexEdge(v, k, oe);
        const auto l{edgeGidToLid[oe]};
        if(oe >= e || l == -1) {
          continue;
        }
        ttk::UnionFind::makeUnion(&uf[i], &uf[l]);
      }
    }
  }

  size_t nLocConnComps{};
  // find connected component ids
  for(size_t i = 0; i < uf.size(); ++i) {
    const auto root{uf[i].find()};
    if(root != nullptr) {
      const auto rank{root->getRank()};
      if(connCompIds[rank] == -1) {
        connCompIds[rank] = nLocConnComps++;
      }
      connCompIds[i] = connCompIds[rank];
    }
  }

  // store component edges ids (useful to get arcs)
  connCompEdges.resize(nLocConnComps);
  for(auto &vec : connCompEdges) {
    // random guess...
    vec.reserve(bucketEdges.size() / nLocConnComps + 1);
  }

  for(size_t i = 0; i < bucketEdges.size(); ++i) {
    const auto e{bucketEdges[i]};
    const auto compId{connCompIds[i]};
    if(compId == -1) {
      continue;
    }
    connCompEdges[compId].emplace_back(e);
    // also put component id on (relevant) edge vertices
    for(SimplexId j = 0; j < 2; ++j) {
      SimplexId v{};
      triangulation.getEdgeVertex(e, j, v);
      if(vertsBucket[v] == bucketId) {
        connComps[v] = compId;
      }
    }
  }

  // cleanup
  for(const auto e : bucketEdges) {
    edgeGidToLid[e] = -1;
  }
  bucketEdges.clear();

  return 0;
}

template <typename triangulationType>
int ttk::Mapper::computeCompBarycenter(
  std::array<float, 3> &baryCoords,
  std::vector<double> &highDimBaryCoords,
  const std::vector<double> &inputHighDimCoords,
  const std::vector<SimplexId> &connCompEdges,
  const triangulationType &triangulation) const {

  size_t nVertsInComp{};

  const auto nVerts{triangulation.getNumberOfVertices()};
  const auto nHighDims{inputHighDimCoords.size() / nVerts};

  for(const auto e : connCompEdges) {
    for(SimplexId j = 0; j < 2; ++j) {
      SimplexId v{};
      triangulation.getEdgeVertex(e, j, v);
      // get barycenter from all vertices from all considered edges
      std::array<float, 3> pt{};
      triangulation.getVertexPoint(v, pt[0], pt[1], pt[2]);
      baryCoords[0] += pt[0];
      baryCoords[1] += pt[1];
      baryCoords[2] += pt[2];
      nVertsInComp++;

      if(this->ComputeHighDimBarycenters) {
        highDimBaryCoords.resize(nHighDims);
        for(size_t i = 0; i < nHighDims; ++i) {
          highDimBaryCoords[i] += inputHighDimCoords[v * nHighDims + i];
        }
      }
    }
  }

  // store barycenter coordinates
  baryCoords[0] /= nVertsInComp;
  baryCoords[1] /= nVertsInComp;
  baryCoords[2] /= nVertsInComp;

  if(this->ComputeHighDimBarycenters) {
    for(size_t i = 0; i < nHighDims; ++i) {
      highDimBaryCoords[i] /= nVertsInComp;
    }
  }

  return 0;
}
