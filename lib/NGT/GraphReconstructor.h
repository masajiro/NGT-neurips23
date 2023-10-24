//
// Copyright (C) 2015 Yahoo Japan Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#pragma once

#include	<unordered_map>
#include	<unordered_set>
#include	<list>

#ifdef _OPENMP
#include	<omp.h>
#else
#warning "*** OMP is *NOT* available! ***"
#endif

#define NGT_SHORTCUT_REDUCTION_WITH_ANGLE
#define NGT_SHORTCUT_REDUCTION_WITH_ADDITIONAL_CONDITION

namespace NGT {

class GraphReconstructor {
 public:
  static void extractGraph(std::vector<NGT::ObjectDistances> &graph, NGT::GraphIndex &graphIndex) {
    graph.reserve(graphIndex.repository.size());
    for (size_t id = 1; id < graphIndex.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "GraphReconstructor::extractGraph: Processed " << id << " objects." << std::endl;
      }
      try {
	NGT::GraphNode &node = *graphIndex.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	NGT::ObjectDistances nd;
	nd.reserve(node.size());
	for (auto n = node.begin(graphIndex.repository.allocator); n != node.end(graphIndex.repository.allocator); ++n) {
	  nd.push_back(ObjectDistance((*n).id, (*n).distance));
        }
	graph.push_back(nd);
#else
	graph.push_back(node);
#endif
	if (graph.back().size() != graph.back().capacity()) {
	  std::cerr << "GraphReconstructor::extractGraph: Warning! The graph size must be the same as the capacity. " << id << std::endl;
	}
      } catch(NGT::Exception &err) {
	graph.push_back(NGT::ObjectDistances());
	continue;
      }
    }

  }




  static void
    adjustPaths(NGT::Index &outIndex)
  {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    std::cerr << "construct index is not implemented." << std::endl;
    exit(1);
#else
    NGT::GraphIndex	&outGraph = dynamic_cast<NGT::GraphIndex&>(outIndex.getIndex());
    size_t rStartRank = 0;
    std::list<std::pair<size_t, NGT::GraphNode> > tmpGraph;
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      NGT::GraphNode &node = *outGraph.getNode(id);
      tmpGraph.push_back(std::pair<size_t, NGT::GraphNode>(id, node));
      if (node.size() > rStartRank) {
	node.resize(rStartRank);
      }
    }
    size_t removeCount = 0;
    for (size_t rank = rStartRank; ; rank++) {
      bool edge = false;
      Timer timer;
      for (auto it = tmpGraph.begin(); it != tmpGraph.end();) {
	size_t id = (*it).first;
	try {
	  NGT::GraphNode &node = (*it).second;
	  if (rank >= node.size()) {
	    it = tmpGraph.erase(it);
	    continue;
	  }
	  edge = true;
	  if (rank >= 1 && node[rank - 1].distance > node[rank].distance) {
	    std::cerr << "distance order is wrong!" << std::endl;
	    std::cerr << id << ":" << rank << ":" << node[rank - 1].id << ":" << node[rank].id << std::endl;
	  }
	  NGT::GraphNode &tn = *outGraph.getNode(id);
	  volatile bool found = false;
	  if (rank < 1000) {
	    for (size_t tni = 0; tni < tn.size() && !found; tni++) {
	      if (tn[tni].id == node[rank].id) {
		continue;
	      }
	      NGT::GraphNode &dstNode = *outGraph.getNode(tn[tni].id);
	      for (size_t dni = 0; dni < dstNode.size(); dni++) {
		if ((dstNode[dni].id == node[rank].id) && (dstNode[dni].distance < node[rank].distance)) {
		  found = true;
		  break;
		}
	      }
	    }
	  } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(10)
#endif
	    for (size_t tni = 0; tni < tn.size(); tni++) {
	      if (found) {
		continue;
	      }
	      if (tn[tni].id == node[rank].id) {
		continue;
	      }
	      NGT::GraphNode &dstNode = *outGraph.getNode(tn[tni].id);
	      for (size_t dni = 0; dni < dstNode.size(); dni++) {
		if ((dstNode[dni].id == node[rank].id) && (dstNode[dni].distance < node[rank].distance)) {
		  found = true;
		}
	      }
	    }
	  }
	  if (!found) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    outGraph.addEdge(id, node.at(i, outGraph.repository.allocator).id,
			     node.at(i, outGraph.repository.allocator).distance, true);
#else
	    tn.push_back(NGT::ObjectDistance(node[rank].id, node[rank].distance));
#endif
	  } else {
	    removeCount++;
	  }
	} catch(NGT::Exception &err) {
	  std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	  it++;
	  continue;
	}
	it++;
      }
      if (edge == false) {
	break;
      }
    }
#endif // NGT_SHARED_MEMORY_ALLOCATOR
  }

  static void
    adjustPathsEffectively(NGT::Index &outIndex, size_t minNoOfEdges = 0)
  {
    NGT::GraphIndex	&outGraph = dynamic_cast<NGT::GraphIndex&>(outIndex.getIndex());
    adjustPathsEffectively(outGraph, minNoOfEdges);
  }

  static bool edgeComp(NGT::ObjectDistance a, NGT::ObjectDistance b) {
    return (a.id & 0x7FFFFFFF) < (b.id & 0x7FFFFFFF);
  }

#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
  static void insert(NGT::GraphNode &node, size_t edgeID, NGT::Distance edgeDistance, NGT::GraphIndex &graph) {
    NGT::ObjectDistance edge(edgeID, edgeDistance);
    GraphNode::iterator ni = std::lower_bound(node.begin(graph.repository.allocator), node.end(graph.repository.allocator), edge, edgeComp);
    node.insert(ni, edge, graph.repository.allocator);
  }

  static bool hasEdge(NGT::GraphIndex &graph, size_t srcNodeID, size_t dstNodeID)
  {
     NGT::GraphNode &srcNode = *graph.getNode(srcNodeID);
     GraphNode::iterator ni = std::lower_bound(srcNode.begin(graph.repository.allocator), srcNode.end(graph.repository.allocator), ObjectDistance(dstNodeID, 0.0), edgeComp);
     return (ni != srcNode.end(graph.repository.allocator)) && ((*ni).id == dstNodeID);
  }

  static bool hasEdge2(NGT::GraphIndex &graph, size_t srcNodeID, size_t dstNodeID, size_t rank = std::numeric_limits<size_t>::max())
  {
     NGT::GraphNode &srcNode = *graph.getNode(srcNodeID);
     auto size = std::min(srcNode.size(), rank + 1);
     for (size_t idx = 0; idx < size; idx++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
       if (srcNode.at(idx, graph.repository.allocator).id == dstNodeID) {
#else
       if (srcNode[idx].id == dstNodeID) {
#endif
	 return true;
       }
     }
     return false;
  }
#else
  static void insert(NGT::GraphNode &node, size_t edgeID, NGT::Distance edgeDistance) {
    NGT::ObjectDistance edge(edgeID, edgeDistance);
    GraphNode::iterator ni = std::lower_bound(node.begin(), node.end(), edge, edgeComp);
    node.insert(ni, edge);
  }

  static bool hasEdge(NGT::GraphIndex &graph, size_t srcNodeID, size_t dstNodeID)
  {
     NGT::GraphNode &srcNode = *graph.getNode(srcNodeID);
     GraphNode::iterator ni = std::lower_bound(srcNode.begin(), srcNode.end(), ObjectDistance(dstNodeID, 0.0), edgeComp);
     return (ni != srcNode.end()) && ((*ni).id == dstNodeID);
  }

  static bool hasEdge2(NGT::GraphIndex &graph, size_t srcNodeID, size_t dstNodeID, size_t rank = std::numeric_limits<size_t>::max())
  {
     NGT::GraphNode &srcNode = *graph.getNode(srcNodeID);
     auto size = std::min(srcNode.size(), rank + 1);
     for (size_t idx = 0; idx < size; idx++) {
       if (srcNode[idx].id == dstNodeID) {
	 return true;
       }
     }
     return false;
  }
#endif


  static void
    adjustPathsEffectively(NGT::GraphIndex &outGraph,
			   size_t minNoOfEdges)
  {
    Timer timer;
    timer.start();
    std::vector<NGT::GraphNode> tmpGraph;
    tmpGraph.reserve(outGraph.repository.size());
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "GraphReconstructor::adjustPaths: # of the extracted nodes=" << id << " peak vm size=" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
	tmpGraph.push_back(node);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	node.clear(outGraph.repository.allocator);
#else
	node.clear();
#endif
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	tmpGraph.push_back(NGT::GraphNode(outGraph.repository.allocator));
#else
	tmpGraph.push_back(NGT::GraphNode());
#endif
      }
    }
    if (outGraph.repository.size() != tmpGraph.size() + 1) {
      std::stringstream msg;
      msg << "GraphReconstructor: Fatal inner error. " << outGraph.repository.size() << ":" << tmpGraph.size() << ", " << outGraph.getPath();
      NGTThrowException(msg);
    }
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: graph preparing time=" << timer << std::endl;
    timer.reset();
    timer.start();

    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> removeCandidates(tmpGraph.size());
    int removeCandidateCount = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t idx = 0; idx < tmpGraph.size(); ++idx) {
      auto it = tmpGraph.begin() + idx;
      size_t id = idx + 1;
      try {
	NGT::GraphNode &srcNode = *it;	
	std::unordered_map<uint32_t, std::pair<uint32_t, float>> neighbors;
	for (uint32_t sni = 0; sni < srcNode.size(); ++sni) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  neighbors[srcNode.at(sni, outGraph.repository.allocator).id] = std::pair<uint32_t, float>(sni, srcNode.at(sni, outGraph.repository.allocator).distance);
#else
	  neighbors[srcNode[sni].id] = std::pair<uint32_t, float>(sni, srcNode[sni].distance);
#endif
	}

	std::vector<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>> candidates;
	for (size_t sni = 0; sni < srcNode.size(); sni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  NGT::GraphNode &pathNode = tmpGraph[srcNode.at(sni, outGraph.repository.allocator).id - 1];
#else
	  NGT::GraphNode &pathNode = tmpGraph[srcNode[sni].id - 1];
#endif
	  for (size_t pni = 0; pni < pathNode.size(); pni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    auto dstNodeID = pathNode.at(pni, outGraph.repository.allocator).id;
#else
	    auto dstNodeID = pathNode[pni].id;
#endif
	    auto dstNode = neighbors.find(dstNodeID);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    if (dstNode != neighbors.end()
		&& srcNode.at(sni, outGraph.repository.allocator).distance < (*dstNode).second.second
		&& pathNode.at(pni, outGraph.repository.allocator).distance < (*dstNode).second.second
		) {
#else
	    if (dstNode != neighbors.end()
		&& srcNode[sni].distance < (*dstNode).second.second
		&& pathNode[pni].distance < (*dstNode).second.second
		) {
#endif
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      candidates.push_back(std::pair<uint32_t, std::pair<uint32_t, uint32_t> >((*dstNode).second.first, std::pair<uint32_t, uint32_t>(srcNode.at(sni, outGraph.repository.allocator).id, dstNodeID)));
#else
	      candidates.push_back(std::pair<uint32_t, std::pair<uint32_t, uint32_t> >((*dstNode).second.first, std::pair<uint32_t, uint32_t>(srcNode[sni].id, dstNodeID)));
#endif
	      removeCandidateCount++;
	    }
	  }
	}
	sort(candidates.begin(), candidates.end(), std::greater<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>>());
	removeCandidates[id - 1].reserve(candidates.size());
	for (size_t i = 0; i < candidates.size(); i++) {
	  removeCandidates[id - 1].push_back(candidates[i].second);
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: extracting removed edge candidates time=" << timer << std::endl;
    std::cerr << "removeCandidateCount=" << removeCandidateCount << std::endl;
    timer.reset();
    timer.start();

    std::list<uint32_t> ids;
    for (size_t idx = 0; idx < tmpGraph.size(); ++idx) {
      ids.push_back(idx + 1);
    }

    int removeCount = 0;
    removeCandidateCount = 0;
    for (size_t rank = 0; ids.size() != 0; rank++) {
      for (auto it = ids.begin(); it != ids.end(); ) {
	size_t id = *it;
	size_t idx = id - 1;
	try {
	  NGT::GraphNode &srcNode = tmpGraph[idx];
	  if (rank >= srcNode.size()) {
	    if (!removeCandidates[idx].empty() && minNoOfEdges == 0) {
	      std::cerr << "Something wrong! ID=" << id << " # of remaining candidates=" << removeCandidates[idx].size() << std::endl;
	      abort();
	    }
#if !defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    NGT::GraphNode empty;
            tmpGraph[idx] = empty;
#endif
	    it = ids.erase(it);
	    continue;
	  }
	  if (removeCandidates[idx].size() > 0 && ((*outGraph.getNode(id)).size() + srcNode.size() - rank) > minNoOfEdges) {
	    removeCandidateCount++;
	    bool pathExist = false;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
            while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == srcNode.at(rank, outGraph.repository.allocator).id)) {
#else
	    while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == srcNode[rank].id)) {
#endif
	      size_t path = removeCandidates[idx].back().first;
	      size_t dst = removeCandidates[idx].back().second;
	      removeCandidates[idx].pop_back();
 	      if (removeCandidates[idx].empty()) {
 	        std::vector<std::pair<uint32_t, uint32_t>> empty;
 		removeCandidates[idx] = empty;
 	      }
              if ((hasEdge(outGraph, id, path)) && (hasEdge(outGraph, path, dst))) {
		pathExist = true;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	        while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == srcNode.at(rank, outGraph.repository.allocator).id)) {
#else
		while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == srcNode[rank].id)) {
#endif
	          removeCandidates[idx].pop_back();
 	          if (removeCandidates[idx].empty()) {
 	            std::vector<std::pair<uint32_t, uint32_t>> empty;
 		    removeCandidates[idx] = empty;
 	          }
		}
		break;
	      }
	    }
	    if (pathExist) {
	      removeCount++;
              it++;
	      continue;
	    }
	  }
	  NGT::GraphNode &outSrcNode = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  insert(outSrcNode, srcNode.at(rank, outGraph.repository.allocator).id, srcNode.at(rank, outGraph.repository.allocator).distance, outGraph);
#else
	  insert(outSrcNode, srcNode[rank].id, srcNode[rank].distance);
#endif
	} catch(NGT::Exception &err) {
	  std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
          it++;
	  continue;
	}
        it++;
      }
    }
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	std::sort(node.begin(outGraph.repository.allocator), node.end(outGraph.repository.allocator));
#else
	std::sort(node.begin(), node.end());
#endif
      } catch(...) {}
    }
  }

  static void removeShortcutEdgesInPartialRanksForTwoStageGraphs(NGT::GraphIndex &outGraph, std::vector<NGT::GraphNode> &tmpGraph,
								 size_t minNoOfEdges, uint32_t rankBegin, uint32_t rankEnd,
								 std::list<uint32_t> &ids, std::vector<bool> &finishedNodes,
								 int &removeCount, int &removeCandidateCount)
  {
    std::cerr << "adjustPathInPartially: rank=" << rankBegin << "->" << rankEnd << std::endl;
    std::cerr << "  vm size(1)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    Timer timer;
    timer.start();

    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> removeCandidates(tmpGraph.size());
    
    std::cerr << "  vm size(2)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    std::cerr << "# of threads=" << omp_get_max_threads() << std::endl;
    auto nthreads = omp_get_max_threads();
    size_t counter[nthreads];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t idx = 0; idx < tmpGraph.size(); ++idx) {
      auto thdID = omp_get_thread_num();
      counter[thdID]++;
      if (counter[thdID] % 100000 == 0) {
	size_t tcount = 0;
	for (int nt = 0; nt != nthreads; nt++) {
	  tcount += counter[nt];
	}
	std::cerr << tcount << " objects have been processed." << std::endl
		  << " vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
      if (finishedNodes[idx]) {
	continue;
      }
      size_t id = idx + 1;
      auto git = tmpGraph.begin() + idx;

      if (id % 100000 == 0) {
	timer.stop();
	std::cerr << "processed " << id << " nodes. # of removed candidates=" << removeCandidateCount << " time=" << timer
		  << " vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
	timer.restart();
      }

      try {
	NGT::GraphNode &srcNode = *git;	
	std::unordered_map<uint32_t, std::pair<uint32_t, float>> dstNodes;
	for (size_t sni = rankBegin; (sni < rankEnd) && (sni < srcNode.size()); sni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  dstNodes[srcNode.at(sni, outGraph.repository.allocator).id] = std::pair<uint32_t, float>(sni, srcNode.at(sni, outGraph.repository.allocator).distance);
#else
	  dstNodes[srcNode[sni].id] = std::pair<uint32_t, float>(sni, srcNode[sni].distance);
#endif
	}
	std::vector<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>> candidates;
	for (size_t sni = 0; sni < srcNode.size(); sni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  auto pathNodeID = srcNode.at(sni, outGraph.repository.allocator).id;
#else
	  auto pathNodeID = srcNode[sni].id;
#endif
	  NGT::GraphNode &pathNode = tmpGraph[pathNodeID - 1];
	  for (size_t pni = 0; pni < pathNode.size(); pni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    auto dstNodeID = pathNode.at(pni, outGraph.repository.allocator).id;
#else
	    auto dstNodeID = pathNode[pni].id;
#endif
	    auto dstNode = dstNodes.find(dstNodeID);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    if (dstNode != dstNodes.end()
		//&& (*pathNode).second.second < srcNode.at(sni, outGraph.repository.allocator).distance
		//&& pathNode.at(pni, outGraph.repository.allocator).distance < srcNode.at(sni, outGraph.repository.allocator).distance
		&& srcNode.at(sni, outGraph.repository.allocator).distance < (*dstNode).second.second
		&& pathNode.at(pni, outGraph.repository.allocator).distance < (*dstNode).second.second
		) {
#else
	    if (dstNode != dstNodes.end()
		//&& (*pathNode).second.second < srcNode[sni].distance
		//&& pathNode[pni].distance < srcNode[sni].distance
		&& srcNode[sni].distance < (*dstNode).second.second
		&& pathNode[pni].distance < (*dstNode).second.second
		) {
#endif
	      candidates.push_back(std::pair<uint32_t, std::pair<uint32_t, uint32_t>>((*dstNode).second.first, std::pair<uint32_t, uint32_t>(pathNodeID, dstNodeID)));
	    }
	  }
	}
	sort(candidates.begin(), candidates.end(), std::greater<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>>());
	removeCandidates[id - 1].reserve(candidates.size());
	for (size_t i = 0; i < candidates.size(); i++) {
	  removeCandidates[id - 1].push_back(candidates[i].second);
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }

    {
      for (size_t idx = 0; idx < tmpGraph.size(); ++idx) {
	removeCandidateCount += removeCandidates[idx].size();
      }
    }
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: extracting removed edge candidates time=" << timer << std::endl;
    std::cerr << "  vm size(3)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    timer.reset();
    timer.start();

    size_t nOfCheckedNodes = 0;
    for (size_t rank = rankBegin; (rank < rankEnd) && (ids.size() != 0); rank++) {

      if (rank % 10 == 0) {
	timer.stop();
        std::cerr << "rank=" << rank << " " << "The number of removed edges=" << removeCount << ":" << removeCandidateCount << " time=" << timer
		  << " vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
	timer.restart();
      }

      for (auto it = ids.begin(); it != ids.end(); ) {
	size_t id = *it;
	size_t idx = id - 1;
	try {
	  NGT::GraphNode &srcNode = tmpGraph[idx];
	  if (rank >= srcNode.size()) {
	    if (!removeCandidates[idx].empty() && minNoOfEdges == 0) {
	      std::stringstream msg;
	      msg << "Something wrong! ID=" << id << " # of remaining candidates=" << removeCandidates[idx].size();
	      NGTThrowException(msg);
	    }
	    it = ids.erase(it);
	    finishedNodes[idx] = true;
	    continue;
	  }
	  if (removeCandidates[idx].size() > 0 && ((*outGraph.getNode(id)).size() + srcNode.size() - rank) > minNoOfEdges) {
	    nOfCheckedNodes++;
	    bool pathExist = false;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
            while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == srcNode.at(rank, outGraph.repository.allocator).id)) {
#else
	    while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == srcNode[rank].id)) {
#endif
	      size_t path = removeCandidates[idx].back().first;
	      size_t dst = removeCandidates[idx].back().second;
	      removeCandidates[idx].pop_back();
 	      if (removeCandidates[idx].empty()) {
 	        std::vector<std::pair<uint32_t, uint32_t>> empty;
 		removeCandidates[idx] = empty;
 	      }
              if ((hasEdge(outGraph, id, path)) && (hasEdge(outGraph, path, dst))) {
		pathExist = true;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	        while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == srcNode.at(rank, outGraph.repository.allocator).id)) {
#else
		while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == srcNode[rank].id)) {
#endif
	          removeCandidates[idx].pop_back();
 	          if (removeCandidates[idx].empty()) {
 	            std::vector<std::pair<uint32_t, uint32_t>> empty;
 		    removeCandidates[idx] = empty;
 	          }
		}
		break;
	      }
	    }
	    if (pathExist) {
	      removeCount++;
              it++;
	      continue;
	    }
	  }
	  NGT::GraphNode &outSrcNode = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  insert(outSrcNode, srcNode.at(rank, outGraph.repository.allocator).id, srcNode.at(rank, outGraph.repository.allocator).distance, outGraph);
#else
	  insert(outSrcNode, srcNode[rank].id, srcNode[rank].distance);
#endif
	} catch(NGT::Exception &err) {
	  std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
          it++;
	  continue;
	}
        it++;
      }
    }
    std::cerr << "removeCandidateCount(2)=" << removeCount << "/" << removeCandidateCount << ":" << nOfCheckedNodes++ << std::endl;
    std::cerr << "  vm size(4)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
  }

  static void
    removeShortcutEdgesInPartialRanksForTwoStageGraphs(NGT::GraphIndex &outGraph,
						       size_t minNoOfEdges)
  {
    Timer timer;
    timer.start();
    std::vector<NGT::GraphNode> tmpGraph;
    tmpGraph.reserve(outGraph.repository.size());
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "GraphReconstructor::adjustPathsInStages: # of the extracted nodes=" << id
		  << " vm size=" << NGT::Common::getProcessVmSizeStr()
		  << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
	tmpGraph.push_back(std::move(node));
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	node.clear(outGraph.repository.allocator);
#else
	node.clear();
#endif
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	tmpGraph.push_back(NGT::GraphNode(outGraph.repository.allocator));
#else
	tmpGraph.push_back(NGT::GraphNode());
#endif
      }
    }
    if (outGraph.repository.size() != tmpGraph.size() + 1) {
      std::stringstream msg;
      msg << "GraphReconstructor: Fatal inner error. " << outGraph.repository.size() << ":" << tmpGraph.size();
      NGTThrowException(msg);
    }
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: graph preparing time=" << timer << std::endl;
    timer.reset();
    timer.start();


    std::list<uint32_t> ids;
    for (size_t idx = 0; idx < tmpGraph.size(); ++idx) {
      ids.push_back(idx + 1);
    }
    std::vector<bool> finishedNodes(tmpGraph.size(), false);
    int removeCount = 0;
    int removeCandidateCount = 0;
    size_t rankStep = 1;
    for (size_t rankBegin = 0; ids.size() != 0;) {
      size_t rankEnd = rankBegin + rankStep;
      int rmcnt = 0;
      removeShortcutEdgesInPartialRanksForTwoStageGraphs(outGraph, tmpGraph, minNoOfEdges, rankBegin, rankEnd, ids, finishedNodes, rmcnt, removeCandidateCount);
      removeCount += rmcnt;
      std::cerr << "removeCount=" << rmcnt << "/" << removeCount << std::endl;
      rankBegin = rankEnd;
    }

    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	std::sort(node.begin(outGraph.repository.allocator), node.end(outGraph.repository.allocator));
#else
	std::sort(node.begin(), node.end());
#endif
      } catch(...) {}
    }
  }

  static void adjustPathsPartially2(NGT::GraphIndex &outGraph, 
				   size_t minNoOfEdges, uint32_t rankBegin, uint32_t rankEnd,
				   std::list<uint32_t> &ids, std::vector<bool> &finishedNodes,
				    int &removeCount, int &totalRemoveCandidateCount, int &removeCandidateCount)
  {
    std::cerr << "adjustPathInPartially2: rank=" << rankBegin << "->" << rankEnd << std::endl;
    std::cerr << "  vm size(1)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    Timer timer;
    timer.start();

    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> removeCandidates(outGraph.repository.size());
    
    std::cerr << "  vm size(2)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    std::cerr << "# of threads=" << omp_get_max_threads() << std::endl;
    auto nthreads = omp_get_max_threads();
    std::vector<size_t> counter(nthreads);
    std::vector<size_t> rmCounter(nthreads);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t id = 1; id < outGraph.repository.size(); ++id) {
      auto thdID = omp_get_thread_num();
      counter[thdID]++;
      if (counter[thdID] % 100000 == 0) {
	size_t tcount = 0;
	for (int nt = 0; nt != nthreads; nt++) {
	  tcount += counter[nt];
	}
	std::cerr << tcount << " objects have been processed." << std::endl
		  << " vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
      if (finishedNodes[id]) {
	continue;
      }

      if (id % 100000 == 0) {
	timer.stop();
	std::cerr << "processed " << id << " nodes. # of removed candidates=" << totalRemoveCandidateCount << " time=" << timer
		  << " vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
	timer.restart();
      }

      try {
	//NGT::GraphNode &srcNode = *git;	
	NGT::GraphNode &srcNode = *outGraph.getNode(id);
	std::unordered_map<uint32_t, std::pair<uint32_t, float>> dstNodes;
	for (size_t sni = rankBegin; (sni < rankEnd) && (sni < srcNode.size()); sni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  dstNodes[srcNode.at(sni, outGraph.repository.allocator).id & 0x7FFFFFFF] = std::pair<uint32_t, float>(sni, srcNode.at(sni, outGraph.repository.allocator).distance);
#else
	  dstNodes[srcNode[sni].id & 0x7FFFFFFF] = std::pair<uint32_t, float>(sni, srcNode[sni].distance);
#endif
	}
	std::vector<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>> candidates;
	for (size_t sni = 0; sni < srcNode.size(); sni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  auto pathNodeID = srcNode.at(sni, outGraph.repository.allocator).id & 0x7FFFFFFF;
#else
	  auto pathNodeID = srcNode[sni].id & 0x7FFFFFFF;
#endif
	  NGT::GraphNode &pathNode = *outGraph.getNode(pathNodeID);
	  for (size_t pni = 0; pni < pathNode.size(); pni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    auto dstNodeID = pathNode.at(pni, outGraph.repository.allocator).id & 0x7FFFFFFF;
#else
	    auto dstNodeID = pathNode[pni].id & 0x7FFFFFFF;
#endif
	    auto dstNode = dstNodes.find(dstNodeID);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    if (dstNode != dstNodes.end()
		&& srcNode.at(sni, outGraph.repository.allocator).distance < (*dstNode).second.second
		&& pathNode.at(pni, outGraph.repository.allocator).distance < (*dstNode).second.second
		) {
#else
	    if (dstNode != dstNodes.end()
		&& srcNode[sni].distance < (*dstNode).second.second
		&& pathNode[pni].distance < (*dstNode).second.second
		) {
#endif
	      candidates.push_back(std::pair<uint32_t, std::pair<uint32_t, uint32_t>>((*dstNode).second.first, std::pair<uint32_t, uint32_t>(pathNodeID, dstNodeID)));
	    }
	  }
	}
	sort(candidates.begin(), candidates.end(), std::greater<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>>());
	removeCandidates[id - 1].reserve(candidates.size());
	rmCounter[thdID] += candidates.size();
	for (size_t i = 0; i < candidates.size(); i++) {
	  removeCandidates[id - 1].push_back(candidates[i].second);
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }
    {
      totalRemoveCandidateCount = 0;
      for (int nt = 0; nt != nthreads; nt++) {
	totalRemoveCandidateCount += rmCounter[nt];
      }
    }
    {
      auto rank = rankEnd - 1;
      removeCandidateCount = 0;
      for (size_t id = 1; id < outGraph.repository.size(); ++id) {
	NGT::GraphNode &srcNode = *outGraph.getNode(id);
	if (rank >= srcNode.size()) continue;
	for (auto &e : removeCandidates[id - 1]) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  if ((srcNode.at(rank, outGraph.repository.allocator).id & 0x7FFFFFFF) == e.second) removeCandidateCount++;
#else
	  if ((srcNode[rank].id & 0x7FFFFFFF) == e.second) removeCandidateCount++;
#endif
	}
      }
    }
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: extracting removed edge candidates time=" << timer << std::endl;
    std::cerr << "  vm size(3)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    timer.reset();
    timer.start();

    size_t nOfCheckedNodes = 0;
    for (size_t rank = rankBegin; (rank < rankEnd) && (ids.size() != 0); rank++) {

      if (rank % 10 == 0) {
	timer.stop();
        std::cerr << "rank=" << rank << " " << "The number of removed edges=" << removeCount << ":" << totalRemoveCandidateCount << " time=" << timer
		  << " vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
	timer.restart();
      }

      for (auto it = ids.begin(); it != ids.end(); ) {
	size_t id = *it;
	size_t idx = id - 1;
	try {
	  //NGT::GraphNode &srcNode = tmpGraph[idx];
	  NGT::GraphNode &srcNode = *outGraph.getNode(id);
	  if (rank >= srcNode.size()) {
	    if (!removeCandidates[idx].empty() && minNoOfEdges == 0) {
	      std::stringstream msg;
	      msg << "Something wrong! ID=" << id << " # of remaining candidates=" << removeCandidates[idx].size();
	      NGTThrowException(msg);
	    }
	    it = ids.erase(it);
	    finishedNodes[id] = true;
	    continue;
	  }
	  if (removeCandidates[idx].size() > 0) {
	    nOfCheckedNodes++;
	    bool pathExist = false;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
            while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == (srcNode.at(rank, outGraph.repository.allocator).id & 0x7FFFFFFF))) {
#else
	    while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == (srcNode[rank].id & 0x7FFFFFFF))) {
#endif
	      size_t path = removeCandidates[idx].back().first;
	      size_t dst = removeCandidates[idx].back().second;
	      removeCandidates[idx].pop_back();
 	      if (removeCandidates[idx].empty()) {
 	        std::vector<std::pair<uint32_t, uint32_t>> empty;
 		removeCandidates[idx] = empty;
 	      }
              if ((hasEdge2(outGraph, id, path, rank)) && (hasEdge2(outGraph, path, dst, rank))) {
		pathExist = true;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	        while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == (srcNode.at(rank, outGraph.repository.allocator).id & 0x7FFFFFFF))) {
#else
		while (!removeCandidates[idx].empty() && (removeCandidates[idx].back().second == (srcNode[rank].id & 0x7FFFFFFF))) {
#endif
	          removeCandidates[idx].pop_back();
 	          if (removeCandidates[idx].empty()) {
 	            std::vector<std::pair<uint32_t, uint32_t>> empty;
 		    removeCandidates[idx] = empty;
 	          }
		}
		break;
	      }
	    }
	    if (pathExist) {
	      removeCount++;
              it++;
	      continue;
	    }
	  }
#if 0
	  NGT::GraphNode &outSrcNode = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  insert(outSrcNode, srcNode.at(rank, outGraph.repository.allocator).id, srcNode.at(rank, outGraph.repository.allocator).distance, outGraph);
#else
	  insert(outSrcNode, srcNode[rank].id, srcNode[rank].distance);
#endif
#endif
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  if ((srcNode.at(rank, outGraph.repository.allocator).id & 0x80000000) != 0) {
	    srcNode.at(rank, outGraph.repository.allocator).id &= 0x7FFFFFFF;
#else
	  if ((srcNode[rank].id & 0x80000000) != 0) {
	    srcNode[rank].id &= 0x7FFFFFFF;
#endif
	  } else {
	    std::cerr << "Warning: something wrong to activate. " << std::endl;
	  }
	} catch(NGT::Exception &err) {
	  std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
          it++;
	  continue;
	}
        it++;
      }
    }
    std::cerr << "removeCandidateCount(2)=" << removeCount << "/" << totalRemoveCandidateCount << ":" << nOfCheckedNodes << std::endl;
    std::cerr << "  vm size(4)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
  }

  static void
    adjustPathsInStages2(NGT::GraphIndex &outGraph,
			 size_t minNoOfEdges)
  {
    Timer timer;
    timer.start();
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "GraphReconstructor::adjustPathsInStages: # of the extracted nodes=" << id
		  << " vm size=" << NGT::Common::getProcessVmSizeStr()
		  << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
	for (size_t rank = 1; rank < node.size(); rank++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  auto &edge = node.at(rank, outGraph.repository.allocator);
#else
	  auto &edge = node[rank];
#endif
	  if ((edge.id & 0x80000000) != 0) {
	    std::stringstream msg;
	    msg << "ID is too large to reduce edges. ID=" << edge.id;
	    NGTThrowException(msg);
	  }
	  edge.id |= 0x80000000;
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      }
    }

    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: graph preparing time=" << timer << std::endl;
    timer.reset();
    timer.start();


    std::list<uint32_t> ids;
    for (size_t id = 1; id < outGraph.repository.size(); ++id) {
      ids.push_back(id);
    }
    std::vector<bool> finishedNodes(outGraph.repository.size(), false);
    int removeCount = 0;
    int removeCandidateCount = 0;
    int totalRemoveCandidateCount = 0;
    size_t rankStep = 5;	
    size_t memKSize = NGT::Common::getSystemHWM();
    size_t vmSize = NGT::Common::getProcessVmSize();
    size_t prev = 0;
    float rcSize = 0.0;
    for (size_t rankBegin = 1; ids.size() != 0;) {
      size_t rankEnd = rankBegin + rankStep;
      int rmcnt = 0;
      std::cerr << "Info vm size=" << NGT::Common::getProcessVmSizeStr()
		<< ":" << NGT::Common::getProcessVmPeakStr()
		<< ":" << NGT::Common::getSystemHWMStr() << std::endl;
      adjustPathsPartially2(outGraph, minNoOfEdges, rankBegin, rankEnd, ids, finishedNodes, rmcnt, totalRemoveCandidateCount, removeCandidateCount);
      size_t vmPeakSize = NGT::Common::getProcessVmPeak();
      if (rcSize == 0.0) {
	rcSize = static_cast<float>(vmPeakSize - vmSize) / static_cast<float>(totalRemoveCandidateCount);
      }
      std::cerr << "Info rcSize=" << rcSize << std::endl;
      std::cerr << "Info vm size=" << NGT::Common::getProcessVmSizeStr()
		<< ":" << NGT::Common::getProcessVmPeakStr()
		<< ":" << NGT::Common::getSystemHWMStr() << std::endl;
      removeCount += rmcnt;
      std::cerr << "removeCount=" << rmcnt << "/" << removeCandidateCount << "/"
		<< totalRemoveCandidateCount << "/" << removeCount << std::endl;
      if (removeCandidateCount != 0) {
	size_t requiredSizePerOneRank = removeCandidateCount * rcSize;
	if (prev < requiredSizePerOneRank) {
	  std::cerr << "Info 2" << std::endl;
	  rankStep = (memKSize - vmSize) * 8 / 10 / (requiredSizePerOneRank * 3 / 2);
	} else {
	  std::cerr << "Info no" << std::endl;
	  rankStep = (memKSize - vmSize) * 8 / 10 / requiredSizePerOneRank;
	}
	prev = requiredSizePerOneRank;
      } else {
	std::cerr << "removeCandidateCount is zero. # of ids=" << ids.size() << std::endl;
	rankStep = rankEnd - rankBegin;
      }
      rankStep = rankStep == 0 ? 1 : rankStep;
      std::cerr << "Info rankStep=" << rankStep << std::endl;
      rankBegin = rankEnd;
    }

    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	std::cerr << "Not implemented yet." << std::endl;
	abort();
#else
	node.erase(std::remove_if(node.begin(), node.end(), [](NGT::ObjectDistance &n){ return (n.id & 0x80000000) != 0; }), node.end());
#endif
      } catch(...) {}
    }
  }

  static void removeShortcutEdgesInPartialRanks1(NGT::GraphIndex *outGraph, 
				     int &totalRemoveCandidateCount,
				     FILE *sortedRemovedCandidatesFile)
  {
    std::cerr << "adjustPathInPartially3" << std::endl;
    std::cerr << "  vm size(1)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    Timer timer;
    timer.start();

    std::cerr << "  vm size(2)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    std::cerr << "# of threads=" << omp_get_max_threads() << std::endl;
    auto nthreads = omp_get_max_threads();
    std::vector<size_t> counter(nthreads);
    std::vector<size_t> rmCounter(nthreads);

    std::vector<FILE*> removeCandidatesFiles(nthreads);
    for (size_t i = 0; i < removeCandidatesFiles.size(); i++) {
      removeCandidatesFiles[i] = tmpfile();
    }
    auto repositorySize = outGraph->repository.size();
    std::cerr << "  vm size(3)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    uint32_t maxRank = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t id = 1; id < repositorySize; ++id) {
      auto thdID = omp_get_thread_num();
      counter[thdID]++;
      if (counter[thdID] % 100000 == 0) {
	size_t tcount = 0;
	for (int nt = 0; nt != nthreads; nt++) {
	  tcount += counter[nt];
	}
	timer.stop();
	std::cerr << tcount << " nodes have been processed. vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << " time=" << timer << std::endl;
	timer.restart();
      }
      //auto git = tmpGraph.begin() + idx;
      try {
	NGT::GraphNode &srcNode = *outGraph->getNode(id);
	if (srcNode.size() > maxRank) maxRank = srcNode.size();
	std::unordered_map<uint32_t, std::pair<uint32_t, float>> dstNodes;
	for (size_t sni = 1; sni < srcNode.size(); sni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  dstNodes[srcNode.at(sni, outGraph->repository.allocator).id] = std::pair<uint32_t, float>(sni, srcNode.at(sni, outGraph->repository.allocator).distance);
#else
	  dstNodes[srcNode[sni].id] = std::pair<uint32_t, float>(sni, srcNode[sni].distance);
#endif
	}
	std::vector<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>> candidates;
	for (size_t sni = 0; sni < srcNode.size(); sni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  auto pathNodeID = srcNode.at(sni, outGraph->repository.allocator).id;
#else
	  auto pathNodeID = srcNode[sni].id;
#endif
	  NGT::GraphNode &pathNode = *outGraph->getNode(pathNodeID);
	  for (size_t pni = 0; pni < pathNode.size(); pni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    auto dstNodeID = pathNode.at(pni, outGraph->repository.allocator).id;
#else
	    auto dstNodeID = pathNode[pni].id;
#endif
	    auto dstNode = dstNodes.find(dstNodeID);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    if (dstNode != dstNodes.end()
		&& srcNode.at(sni, outGraph->repository.allocator).distance < (*dstNode).second.second
		&& pathNode.at(pni, outGraph->repository.allocator).distance < (*dstNode).second.second
		) {
#else
	    if (dstNode != dstNodes.end()
		&& srcNode[sni].distance < (*dstNode).second.second
		&& pathNode[pni].distance < (*dstNode).second.second
		) {
#endif
	      candidates.push_back(std::pair<uint32_t, std::pair<uint32_t, uint32_t>>((*dstNode).second.first, std::pair<uint32_t, uint32_t>(pathNodeID, dstNodeID)));
	    }
	  }
	}
	sort(candidates.begin(), candidates.end(), std::greater<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>>());
	rmCounter[thdID] += candidates.size();
	uint32_t csize = candidates.size();
	fwrite(&id, sizeof(id), 1, removeCandidatesFiles[thdID]);
	fwrite(&csize, sizeof(csize), 1, removeCandidatesFiles[thdID]);
	for (size_t i = 0; i < candidates.size(); i++) {
	  auto ri = candidates.size() - i - 1;
	  uint32_t r = candidates[ri].first;
	  fwrite(&r, sizeof(r), 1, removeCandidatesFiles[thdID]);
	  fwrite(&candidates[ri].second.first, sizeof(candidates[ri].second.first), 1, removeCandidatesFiles[thdID]);
	  fwrite(&candidates[ri].second.second, sizeof(candidates[ri].second.second), 1, removeCandidatesFiles[thdID]);
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }

    std::cerr << "maxRank=" << maxRank << std::endl;
    std::cerr << "  vm size(4)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    delete outGraph;
    outGraph = 0;
    std::cerr << "delete  vm size(5)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
   {
     std::cerr << "new  vm size(6)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
     std::vector<std::pair<size_t, size_t>> pos(repositorySize);
     std::vector<uint32_t> fdTable;
     std::cerr << "vm size(6.1)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
     for (size_t fi = 0; fi < removeCandidatesFiles.size(); fi++) {
       std::cerr << "FI=" << fi << std::endl;
       fseek(removeCandidatesFiles[fi], 0, SEEK_SET);
       uint32_t cid, csize;
       for (uint32_t id = 1; id < repositorySize; ++id) {
	 auto n = fread(&cid, sizeof(cid), 1, removeCandidatesFiles[fi]);
	 if (n == 0) break;
	 fread(&csize, sizeof(csize), 1, removeCandidatesFiles[fi]);
	 if (std::lower_bound(fdTable.begin(), fdTable.end(), cid) != fdTable.end()) {
	   std::cerr << "Fatal error!." << std::endl;
	   abort();
	 }
	 if (fdTable.size() <= fi) {
	   fdTable.push_back(cid);
	 }
	 pos[cid].first = ftell(removeCandidatesFiles[fi]);
	 pos[cid].second = pos[cid].first + csize * sizeof(uint32_t) * 3;
	 fseek(removeCandidatesFiles[fi], pos[cid].second, SEEK_SET);
       }
     }
     std::cerr << "vm size(6.2)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
     std::cerr << "end of initialize." << std::endl;
     uint32_t rankStep = 10;
     for (uint32_t rcRank = 1; rcRank < maxRank; rcRank += rankStep) {
       std::cerr << "rank=" << rcRank << std::endl;
       std::vector<std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>>> rm(repositorySize);
       for (uint32_t id = 1; id < repositorySize; ++id) {
	 auto fd = std::upper_bound(fdTable.begin(), fdTable.end(), id);
	 auto fidx = std::distance(fdTable.begin(), fd) - 1;
	 while (pos[id].first != pos[id].second) {
	   fseek(removeCandidatesFiles[fidx], pos[id].first, SEEK_SET);
	   uint32_t r, path, dst;
	   fread(&r, sizeof(r), 1, removeCandidatesFiles[fidx]);
	   fread(&path, sizeof(path), 1, removeCandidatesFiles[fidx]);
	   fread(&dst, sizeof(dst), 1, removeCandidatesFiles[fidx]);
	   if (r < rcRank + rankStep) {
	     rm[id].push_back(std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>(r, id, path, dst));
	     pos[id].first = ftell(removeCandidatesFiles[fidx]);
	   } else {
	     break;
	   }
	 }
       }
       std::cerr << "end of extracting rm." << std::endl;
       for (uint32_t id = 1; id < repositorySize; ++id) {
	 std::reverse(rm[id].begin(), rm[id].end());
       }
       std::cerr << "end of reverse." << std::endl;
       std::cerr << "vm size(8)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
       for (auto r = rcRank; r < rcRank + rankStep; r++) {
	 for (uint32_t id = 1; id < repositorySize; ++id) {
	   if (rm[id].empty()) continue;
	   if (std::get<1>(rm[id].back()) != id) {
	     std::cerr << "Fatal Error ID. " << id << ":" << std::get<0>(rm[id].back()) << std::endl;
	     abort();
	   }
	   while (std::get<0>(rm[id].back()) == r) {
	     fwrite(&std::get<0>(rm[id].back()), sizeof(uint32_t), 1, sortedRemovedCandidatesFile);
	     fwrite(&std::get<1>(rm[id].back()), sizeof(uint32_t), 1, sortedRemovedCandidatesFile);
	     fwrite(&std::get<2>(rm[id].back()), sizeof(uint32_t), 1, sortedRemovedCandidatesFile);
	     fwrite(&std::get<3>(rm[id].back()), sizeof(uint32_t), 1, sortedRemovedCandidatesFile);
	     rm[id].pop_back();
	   }
	 }
       }
       for (size_t i = 0; i < rm.size(); i++) {
	 if (!rm[i].empty()) {
	   std::cerr << "Not empty. " << rm[i].size() << std::endl;
	 }
       }
     }
   }

    for (auto &removeCandidatesFile : removeCandidatesFiles) {
      fclose(removeCandidatesFile);
    }
    {
      totalRemoveCandidateCount = 0;
      for (int nt = 0; nt != nthreads; nt++) {
	totalRemoveCandidateCount += rmCounter[nt];
      }
    }
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: extracting removed edge candidates time=" << timer << std::endl;
    std::cerr << "  vm size(3)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    std::cerr << "removeCandidateCount(1.5)=" << totalRemoveCandidateCount << std::endl;

    timer.reset();
    timer.start();

  }

  static void removeShortcutEdgesInPartialRanks2(NGT::GraphIndex &outGraph, 
				     std::list<uint32_t> &ids,
				     int &removeCount,
				     FILE *sortedRemovedCandidatesFile)
  {
    Timer timer;

    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "GraphReconstructor::adjustPathsInStages: # of the extracted nodes=" << id
		  << " vm size=" << NGT::Common::getProcessVmSizeStr()
		  << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
	for (size_t rank = 1; rank < node.size(); rank++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  auto &edge = node.at(rank, outGraph.repository.allocator);
#else
	  auto &edge = node[rank];
#endif
	  if ((edge.id & 0x80000000) != 0) {
	    std::stringstream msg;
	    msg << "ID is too large to reduce edges. ID=" << edge.id;
	    NGTThrowException(msg);
	  }
	  edge.id |= 0x80000000;
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      }
    }

    fseek(sortedRemovedCandidatesFile, 0, SEEK_SET);
    uint32_t crank, cid, cpath, cdst;
    fread(&crank, sizeof(crank), 1, sortedRemovedCandidatesFile);
    fread(&cid, sizeof(cid), 1, sortedRemovedCandidatesFile);
    fread(&cpath, sizeof(cpath), 1, sortedRemovedCandidatesFile);
    fread(&cdst, sizeof(cdst), 1, sortedRemovedCandidatesFile);

    size_t nOfCheckedNodes = 0;
    for (size_t rank = 1; ids.size() != 0; rank++) {

      if (rank % 10 == 0) {
	timer.stop(); 
       std::cerr << "rank=" << rank << " " << "The number of removed edges=" << removeCount << " time=" << timer
		  << " vm size=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
	timer.restart();
      }

      for (auto it = ids.begin(); it != ids.end(); ) {
	size_t id = *it;
	try {
	  NGT::GraphNode &srcNode = *outGraph.getNode(id);
	  if (rank >= srcNode.size()) {
	    it = ids.erase(it);
	    continue;
	  }
	  bool pathExist = false;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  while ((rank == crank) && (id == cid) && (cdst == (srcNode.at(rank, outGraph.repository.allocator).id & 0x7FFFFFFF))) {
#else
	  while ((rank == crank) && (id == cid) && (cdst == (srcNode[rank].id & 0x7FFFFFFF))) {

#endif
	    nOfCheckedNodes++;
	    size_t path = cpath;
	    size_t dst = cdst;
	    if (id != cid) {
	      std::cerr << "Wrong id " << id << ":" << cid << std::endl;
	      abort();
	    }
	    if (rank != crank) {
	      std::cerr << "Wrong rank " << rank << ":" << crank << std::endl;
	      abort();
	    }
	    {
	      auto n = fread(&crank, sizeof(crank), 1, sortedRemovedCandidatesFile);
	      if (n == 0) {
		cid = 0;
		std::cerr << "End of file." << std::endl;
	      } else {
		fread(&cid, sizeof(cid), 1, sortedRemovedCandidatesFile);
		fread(&cpath, sizeof(cpath), 1, sortedRemovedCandidatesFile);
		fread(&cdst, sizeof(cdst), 1, sortedRemovedCandidatesFile);
	      }
	    }
	    if ((hasEdge2(outGraph, id, path, rank)) && (hasEdge2(outGraph, path, dst, rank))) {
	      pathExist = true;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      while ((rank == crank) && (id == cid) && (cdst == (srcNode.at(rank, outGraph.repository.allocator).id & 0x7FFFFFFF))) {
#else
	      while ((rank == crank) && (id == cid) && (cdst == (srcNode[rank].id & 0x7FFFFFFF))) {
#endif
		{
		  auto n = fread(&crank, sizeof(crank), 1, sortedRemovedCandidatesFile);
		  if (n == 0) {
		    cid = 0;
		    std::cerr << "End of file." << std::endl;
		    break;
		  } else {
		    fread(&cid, sizeof(cid), 1, sortedRemovedCandidatesFile);
		    fread(&cpath, sizeof(cpath), 1, sortedRemovedCandidatesFile);
		    fread(&cdst, sizeof(cdst), 1, sortedRemovedCandidatesFile);
		  }
		}
	      }
	      break;
	    }
	  }
	  if (pathExist) {
	    removeCount++;
	    it++;
	    continue;
	  }
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  if ((srcNode.at(rank, outGraph.repository.allocator).id & 0x80000000) != 0) {
	    srcNode.at(rank, outGraph.repository.allocator).id &= 0x7FFFFFFF;
#else
	  if ((srcNode[rank].id & 0x80000000) != 0) {
	    srcNode[rank].id &= 0x7FFFFFFF;
#endif
	  } else {
	    std::cerr << "Warning: something wrong to activate. " << std::endl;
	  }
	} catch(NGT::Exception &err) {
	  std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
          it++;
	  continue;
	}
        it++;
      }
    }

    std::cerr << "removeCandidateCount(2)=" << removeCount << ":" << nOfCheckedNodes << std::endl;
    std::cerr << "  vm size(4)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
  }

  static void
    removeShortcutEdgesInPartialRanks(NGT::GraphIndex **outGraph, const std::string &outIndexPath,
				      size_t minNoOfEdges)
  {
    Timer timer;
    timer.start();
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: graph preparing time=" << timer << std::endl;
    timer.reset();
    timer.start();

    int removeCount = 0;
    int removeCandidateCount = 0;

    std::cerr << "Info index path=" << outIndexPath << std::endl;
    std::cerr << "Info vm size=" << NGT::Common::getProcessVmSizeStr()
	      << ":" << NGT::Common::getProcessVmPeakStr()
	      << ":" << NGT::Common::getSystemHWMStr() << std::endl;
    FILE *sortedRemovedCandidatesFile = tmpfile();
    removeShortcutEdgesInPartialRanks1(*outGraph,  removeCandidateCount, sortedRemovedCandidatesFile);

#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    *outGraph = new NGT::GraphIndex(outIndexPath, false);
#else
    *outGraph = new NGT::GraphIndex(outIndexPath, false, NGT::Index::OpenTypeObjectDisabled);
#endif
    std::cerr << "Reopen" << std::endl;
    std::cerr << "GraphOptimizer::execute:  vm size=" << NGT::Common::getProcessVmSizeStr()
	      << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    std::cerr << " delete all of objects" << std::endl;
    try {
      (*outGraph)->destructObjectSpace();
    } catch (NGT::Exception &err) {
      delete *outGraph;
      *outGraph = 0;
      throw(err);
    }
    std::list<uint32_t> ids;
    for (size_t id = 1; id < (*outGraph)->repository.size(); ++id) {
      ids.push_back(id);
    }
    removeShortcutEdgesInPartialRanks2(**outGraph, ids, removeCount, sortedRemovedCandidatesFile);
    std::cerr << "Info vm size=" << NGT::Common::getProcessVmSizeStr()
	      << ":" << NGT::Common::getProcessVmPeakStr()
	      << ":" << NGT::Common::getSystemHWMStr() << std::endl;
    std::cerr << "removeCount=" << removeCandidateCount << "/" << removeCount << std::endl;

    for (size_t id = 1; id < (*outGraph)->repository.size(); id++) {
      try {
	NGT::GraphNode &node = *(*outGraph)->getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	std::cerr << "No implemented yet." << std::endl;
	abort();
#else
	node.erase(std::remove_if(node.begin(), node.end(), [](NGT::ObjectDistance &n){ return (n.id & 0x80000000) != 0; }), node.end());
#endif
      } catch(...) {}
    }
    std::cerr << "end of adjustPathsInStages. the number of remved edges=" << removeCount << std::endl;
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths removing edges time=" << timer << std::endl;
  }

  static void removeShortcutEdges(NGT::GraphIndex &outGraph, int &removeCandidateCount, float range, size_t nOfThreads = 0)
  {
    std::cerr << "removeShortcutEdges";
#ifdef NGT_SHORTCUT_REDUCTION_WITH_ANGLE
    std::cerr << " with angle" << std::endl;
#endif
    std::cerr << ", range=" << range << std::endl;
    std::cerr << "  vm size(1)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    Timer timer;
    timer.start();

    uint32_t maxRank = 0;
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "GraphReconstructor::adjustPathsInStages: # of the extracted nodes=" << id
		  << " vm size=" << NGT::Common::getProcessVmSizeStr()
		  << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
	if (node.size() > maxRank) maxRank = node.size();
	for (size_t rank = 1; rank < node.size(); rank++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  auto &edge = node.at(rank, outGraph.repository.allocator);
#else
	  auto &edge = node[rank];
#endif
	  if ((edge.id & 0x80000000) != 0) {
	    std::stringstream msg;
	    msg << "ID is too large to reduce edges. ID=" << edge.id;
	    NGTThrowException(msg);
	  }
	  edge.id |= 0x80000000;
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
      }
    }

    std::cerr << "  vm size(2)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    std::cerr << "# of max threads=" << omp_get_max_threads() << std::endl;
    if (nOfThreads != 0) {
      omp_set_num_threads(nOfThreads);
    }
    auto nthreads = omp_get_max_threads();
    std::cerr << "# of threads=" << nthreads << std::endl;
    removeCandidateCount = 0;
    auto repositorySize = outGraph.repository.size();
    std::cerr << "  vm size(3)=" << NGT::Common::getProcessVmSizeStr() << ":" << NGT::Common::getProcessVmPeakStr() << std::endl;
    double prevTime = 0.0;
    for (uint32_t rank = 1; rank < maxRank; rank++) {
      timer.stop();
      if (timer.time - prevTime > 4.0) {
	std::cerr << "rank=" << rank << " "
		  << NGT::Common::getProcessVmSizeStr() << ":"
		  << NGT::Common::getProcessVmPeakStr()
		  << " time=" << timer << std::endl;
	prevTime = timer.time;
	timer.restart();
      }
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (uint32_t id = 1; id < repositorySize; ++id) {
	auto thdID = omp_get_thread_num();
	try {
	  NGT::GraphNode &srcNode = *outGraph.getNode(id);
	  if (rank >= srcNode.size()) continue;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  auto dstNodeID = srcNode.at(rank, outGraph.repository.allocator).id & 0x7FFFFFFF;
	  auto dstNodeDistance = srcNode.at(rank, outGraph.repository.allocator).distance;
#else
	  auto dstNodeID = srcNode[rank].id & 0x7FFFFFFF;
	  auto dstNodeDistance = srcNode[rank].distance;
#endif
#ifdef NGT_SHORTCUT_REDUCTION_WITH_ANGLE
	  auto dstNodeDistance2 = dstNodeDistance * dstNodeDistance;
#endif
	  bool found = false;
	  for (size_t sni = 0; sni < srcNode.size() && sni < rank; sni++) {
	  //for (size_t sni = 0; sni < srcNode.size(); sni++) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    auto pathNodeID = srcNode.at(sni, outGraph.repository.allocator).id;
#else
	    auto pathNodeID = srcNode[sni].id;
#endif
	    if ((pathNodeID & 0x80000000) != 0) continue;
#ifdef NGT_SHORTCUT_REDUCTION_WITH_ANGLE
	    auto srcNodeDistance2 = srcNode[sni].distance * srcNode[sni].distance;
#else
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    if (srcNode.at(sni, outGraph.repository.allocator).distance >= dstNodeDistance) continue;
#else
	    if (srcNode[sni].distance >= dstNodeDistance) continue;
#endif
#endif
	    NGT::GraphNode &pathNode = *outGraph.getNode(pathNodeID);
	    for (size_t pni = 0; pni < pathNode.size() && pni <= rank; pni++) {
#ifdef NGT_SHORTCUT_REDUCTION_WITH_ANGLE
	      auto pathNodeDistance2 = pathNode[pni].distance * pathNode[pni].distance;
	      auto v1 = srcNodeDistance2 + pathNodeDistance2 - dstNodeDistance2;
	      auto v2 = 2.0 * srcNode[sni].distance * pathNode[pni].distance;
	      auto cosAlpha = v1 / v2;
	      if (cosAlpha >= range) {
		break;
	      }
#else
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      if (pathNode.at(pni, outGraph.repository.allocator).distance >= dstNodeDistance) break;
#else
	      if (pathNode[pni].distance >= dstNodeDistance) break;
#endif
#ifdef NGT_SHORTCUT_REDUCTION_WITH_ADDITIONAL_CONDITION
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      if (srcNode.at(sni, outGraph.repository.allocator).distance + pathNode.at(pni, outGraph.repository.allocator).distance >= dstNodeDistance * range) break;
#else
	      if (srcNode[sni].distance + pathNode[pni].distance >= dstNodeDistance * range) break;
#endif
#endif
#endif
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      auto nodeID = pathNode.at(pni, outGraph.repository.allocator).id;
#else
	      auto nodeID = pathNode[pni].id;
#endif
	      if ((nodeID & 0x80000000) != 0) continue;
	      if (nodeID != dstNodeID) continue;
	      found = true;
	      removeCandidateCount++;
	      break;
	    }
	    if (found) break;
	  }
	  if (!found) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    srcNode.at(rank, outGraph.repository.allocator).id &= 0x7FFFFFFF;
#else
	    srcNode[rank].id &= 0x7FFFFFFF;
#endif
	  }
	} catch(NGT::Exception &err) {
	  std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	  continue;
	}
      }
    }

    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	std::cerr << "Not implemented yet." << std::endl;
	abort();
#else
	node.erase(std::remove_if(node.begin(), node.end(), [](NGT::ObjectDistance &n){ return (n.id & 0x80000000) != 0; }), node.end());
#endif
      } catch(...) {}
    }

  }

  static void
    removeShortcutEdges(NGT::GraphIndex &outGraph, const std::string &outIndexPath, 
			float range, size_t nOfThreads, size_t minNoOfEdges)
  {
    Timer timer;
    timer.start();
    timer.stop();
    std::cerr << "GraphReconstructor::adjustPaths: graph preparing time=" << timer << std::endl;
    timer.reset();
    timer.start();

    int removeCount = 0;
    int removeCandidateCount = 0;

    std::cerr << "Info vm size=" << NGT::Common::getProcessVmSizeStr()
	      << ":" << NGT::Common::getProcessVmPeakStr()
	      << ":" << NGT::Common::getSystemHWMStr() << std::endl;

    removeShortcutEdges(outGraph, removeCandidateCount, range, nOfThreads);

  }


  static
    void convertToANNG(std::vector<NGT::ObjectDistances> &graph)
  {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    std::cerr << "convertToANNG is not implemented for shared memory." << std::endl;
    return;
#else
    std::cerr << "convertToANNG begin" << std::endl;
    for (size_t idx = 0; idx < graph.size(); idx++) {
      NGT::GraphNode &node = graph[idx];
      for (auto ni = node.begin(); ni != node.end(); ++ni) {
	graph[(*ni).id - 1].push_back(NGT::ObjectDistance(idx + 1, (*ni).distance));
      }
    }
    for (size_t idx = 0; idx < graph.size(); idx++) {
      NGT::GraphNode &node = graph[idx];
      if (node.size() == 0) {
	continue;
      }
      std::sort(node.begin(), node.end());
      NGT::ObjectID prev = 0;
      for (auto it = node.begin(); it != node.end();) {
	if (prev == (*it).id) {
	  it = node.erase(it);
	  continue;
	}
	prev = (*it).id;
	it++;
      }
      NGT::GraphNode tmp = node;
      node.swap(tmp);
    }
    std::cerr << "convertToANNG end" << std::endl;
#endif
  }

  static
    void reconstructGraph(std::vector<NGT::ObjectDistances> &graph, NGT::GraphIndex &outGraph, size_t originalEdgeSize, size_t reverseEdgeSize)
  {
    if (reverseEdgeSize > 10000) {
      std::cerr << "something wrong. Edge size=" << reverseEdgeSize << std::endl;
      exit(1);
    }

    NGT::Timer	originalEdgeTimer, reverseEdgeTimer, normalizeEdgeTimer;
    originalEdgeTimer.start();

    size_t warningCount = 0;
    const size_t warningLimit = 10;
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
	if (originalEdgeSize == 0) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  node.clear(outGraph.repository.allocator);
#else
	  NGT::GraphNode empty;
	  node.swap(empty);
#endif
	} else {
	  NGT::ObjectDistances n = graph[id - 1];
	  if (n.size() < originalEdgeSize) {
	    warningCount++;
	    if (warningCount <= warningLimit) {
	      std::cerr << "GraphReconstructor: Warning. The edges are too few. " << n.size() << ":" << originalEdgeSize << " for " << id << std::endl;
	    }
	    if (warningCount == warningLimit) {
	      std::cerr << "GraphReconstructor: Info. Too many warnings. Warning is disabled." << std::endl;
	    }
	    continue;
	  }
	  n.resize(originalEdgeSize);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	  node.copy(n, outGraph.repository.allocator);
#else
	  node.swap(n);
#endif
	}
      } catch(NGT::Exception &err) {
	warningCount++;
	if (warningCount <= warningLimit) {
	  std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	}
	if (warningCount == warningLimit) {
	  std::cerr << "GraphReconstructor: Info. Too many warnings. Warning is disabled." << std::endl;
	}
	continue;
      }
    }
    if (warningCount > warningLimit) {
      std::cerr << "GraphReconstructor: The total " << warningCount << " Warnings." << std::endl;
    }
    originalEdgeTimer.stop();

    reverseEdgeTimer.start();
    int insufficientNodeCount = 0;
    for (size_t id = 1; id <= graph.size(); ++id) {
      try {
	NGT::ObjectDistances &node = graph[id - 1];
	size_t rsize = reverseEdgeSize;
	if (rsize > node.size()) {
	  insufficientNodeCount++;
	  rsize = node.size();
	}
	for (size_t i = 0; i < rsize; ++i) {
	  NGT::Distance distance = node[i].distance;
	  size_t nodeID = node[i].id;
	  try {
	    NGT::GraphNode &n = *outGraph.getNode(nodeID);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    n.push_back(NGT::ObjectDistance(id, distance), outGraph.repository.allocator);
#else
	    n.push_back(NGT::ObjectDistance(id, distance));
#endif
	  } catch(...) {}
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }
    reverseEdgeTimer.stop();
    if (insufficientNodeCount != 0) {
      std::cerr << "# of the nodes edges of which are in short = " << insufficientNodeCount << std::endl;
    }

    normalizeEdgeTimer.start();
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      try {
	NGT::GraphNode &n = *outGraph.getNode(id);
	if (id % 100000 == 0) {
	  std::cerr << "Processed " << id << " nodes" << std::endl;
	}
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	std::sort(n.begin(outGraph.repository.allocator), n.end(outGraph.repository.allocator));
#else
	std::sort(n.begin(), n.end());
#endif
	NGT::ObjectID prev = 0;
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	for (auto it = n.begin(outGraph.repository.allocator); it != n.end(outGraph.repository.allocator);) {
#else
	for (auto it = n.begin(); it != n.end();) {
#endif
	  if (prev == (*it).id) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    it = n.erase(it, outGraph.repository.allocator);
#else
	    it = n.erase(it);
#endif
	    continue;
	  }
	  prev = (*it).id;
	  it++;
	}
#if !defined(NGT_SHARED_MEMORY_ALLOCATOR)
	NGT::GraphNode tmp = n;
	n.swap(tmp);
#endif
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }
    normalizeEdgeTimer.stop();
    std::cerr << "Reconstruction time=" << originalEdgeTimer.time << ":" << reverseEdgeTimer.time
	 << ":" << normalizeEdgeTimer.time << std::endl;

    NGT::Property prop;
    outGraph.getProperty().get(prop);
    prop.graphType = NGT::NeighborhoodGraph::GraphTypeONNG;
    outGraph.getProperty().set(prop);
  }



  static
    void reconstructGraphWithConstraint(std::vector<NGT::ObjectDistances> &graph, NGT::GraphIndex &outGraph,
					size_t originalEdgeSize, size_t reverseEdgeSize,
					char mode = 'a')
  {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    std::cerr << "reconstructGraphWithConstraint is not implemented." << std::endl;
    abort();
#else

    NGT::Timer	originalEdgeTimer, reverseEdgeTimer, normalizeEdgeTimer;

    if (reverseEdgeSize > 10000) {
      std::cerr << "something wrong. Edge size=" << reverseEdgeSize << std::endl;
      exit(1);
    }

    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "Processed " << id << std::endl;
      }
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
	if (node.size() == 0) {
	  continue;
	}
	node.clear();
	NGT::GraphNode empty;
	node.swap(empty);
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }
    NGT::GraphIndex::showStatisticsOfGraph(outGraph);

    std::vector<ObjectDistances> reverse(graph.size() + 1);	
    for (size_t id = 1; id <= graph.size(); ++id) {
      try {
	NGT::GraphNode &node = graph[id - 1];
	if (id % 100000 == 0) {
	  std::cerr << "Processed (summing up) " << id << std::endl;
	}
	for (size_t rank = 0; rank < node.size(); rank++) {
	  reverse[node[rank].id].push_back(ObjectDistance(id, node[rank].distance));
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }

    std::vector<std::pair<size_t, size_t> > reverseSize(graph.size() + 1);	
    reverseSize[0] = std::pair<size_t, size_t>(0, 0);
    for (size_t rid = 1; rid <= graph.size(); ++rid) {
      reverseSize[rid] = std::pair<size_t, size_t>(reverse[rid].size(), rid);
    }
    std::sort(reverseSize.begin(), reverseSize.end());		


    std::vector<uint32_t> indegreeCount(graph.size(), 0);	
    size_t zeroCount = 0;
    for (size_t sizerank = 0; sizerank <= reverseSize.size(); sizerank++) {

      if (reverseSize[sizerank].first == 0) {
	zeroCount++;
	continue;
      }
      size_t rid = reverseSize[sizerank].second;	
      ObjectDistances &rnode = reverse[rid];		
      for (auto rni = rnode.begin(); rni != rnode.end(); ++rni) {
	if (indegreeCount[(*rni).id] >= reverseEdgeSize) {	
	  continue;
	}
	NGT::GraphNode &node = *outGraph.getNode(rid);	
	if (indegreeCount[(*rni).id] > 0 && node.size() >= originalEdgeSize) {
	  continue;
	}
	
	node.push_back(NGT::ObjectDistance((*rni).id, (*rni).distance));
	indegreeCount[(*rni).id]++;
      }
    }
    reverseEdgeTimer.stop();
    std::cerr << "The number of nodes with zero outdegree by reverse edges=" << zeroCount << std::endl;
    NGT::GraphIndex::showStatisticsOfGraph(outGraph);

    normalizeEdgeTimer.start();
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      try {
	NGT::GraphNode &n = *outGraph.getNode(id);
	if (id % 100000 == 0) {
	  std::cerr << "Processed " << id << std::endl;
	}
	std::sort(n.begin(), n.end());
	NGT::ObjectID prev = 0;
	for (auto it = n.begin(); it != n.end();) {
	  if (prev == (*it).id) {
	    it = n.erase(it);
	    continue;
	  }
	  prev = (*it).id;
	  it++;
	}
	NGT::GraphNode tmp = n;
	n.swap(tmp);
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }
    normalizeEdgeTimer.stop();
    NGT::GraphIndex::showStatisticsOfGraph(outGraph);

    originalEdgeTimer.start();
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "Processed " << id << std::endl;
      }
      NGT::GraphNode &node = graph[id - 1];
      try {
	NGT::GraphNode &onode = *outGraph.getNode(id);
	bool stop = false;
	for (size_t rank = 0; (rank < node.size() && rank < originalEdgeSize) && stop == false; rank++) {
	  switch (mode) {
	  case 'a':
	    if (onode.size() >= originalEdgeSize) {
	      stop = true;
	      continue;
	    }
	    break;
	  case 'c':
	    break;
	  }
	  NGT::Distance distance = node[rank].distance;
	  size_t nodeID = node[rank].id;
	  outGraph.addEdge(id, nodeID, distance, false);
	}
      } catch(NGT::Exception &err) {
	std::cerr << "GraphReconstructor: Warning. Cannot get the node. ID=" << id << ":" << err.what() << std::endl;
	continue;
      }
    }
    originalEdgeTimer.stop();
    NGT::GraphIndex::showStatisticsOfGraph(outGraph);

    std::cerr << "Reconstruction time=" << originalEdgeTimer.time << ":" << reverseEdgeTimer.time
	 << ":" << normalizeEdgeTimer.time << std::endl;

#endif
  }

  // reconstruct a pseudo ANNG with a fewer edges from an actual ANNG with more edges.
  // graph is a source ANNG
  // index is an index with a reconstructed ANNG
  static
    void reconstructANNGFromANNG(std::vector<NGT::ObjectDistances> &graph, NGT::Index &index, size_t edgeSize)
  {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    std::cerr << "reconstructANNGFromANNG is not implemented." << std::endl;
    abort();
#else

    NGT::GraphIndex	&outGraph = dynamic_cast<NGT::GraphIndex&>(index.getIndex());

    // remove all edges in the index.
    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      if (id % 1000000 == 0) {
	std::cerr << "Processed " << id << " nodes." << std::endl;
      }
      try {
	NGT::GraphNode &node = *outGraph.getNode(id);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	node.clear(outGraph.repository.allocator);
#else
	NGT::GraphNode empty;
	node.swap(empty);
#endif
      } catch(NGT::Exception &err) {
      }
    }

    for (size_t id = 1; id <= graph.size(); ++id) {
      size_t edgeCount = 0;
      try {
	NGT::ObjectDistances &node = graph[id - 1];
	NGT::GraphNode &n = *outGraph.getNode(id);
	NGT::Distance prevDistance = 0.0;
	assert(n.size() == 0);
	for (size_t i = 0; i < node.size(); ++i) {
	  NGT::Distance distance = node[i].distance;
	  if (prevDistance > distance) {
	    NGTThrowException("Edge distance order is invalid");
	  }
	  prevDistance = distance;
	  size_t nodeID = node[i].id;
	  if (node[i].id < id) {
	    try {
	      NGT::GraphNode &dn = *outGraph.getNode(nodeID);
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	      n.push_back(NGT::ObjectDistance(nodeID, distance), outGraph.repository.allocator);
	      dn.push_back(NGT::ObjectDistance(id, distance), outGraph.repository.allocator);
#else
	      n.push_back(NGT::ObjectDistance(nodeID, distance));
	      dn.push_back(NGT::ObjectDistance(id, distance));
#endif
	    } catch(...) {}
	    edgeCount++;
	  }
	  if (edgeCount >= edgeSize) {
	    break;
	  }
	}
      } catch(NGT::Exception &err) {
      }
    }

    for (size_t id = 1; id < outGraph.repository.size(); id++) {
      try {
	NGT::GraphNode &n = *outGraph.getNode(id);
	std::sort(n.begin(), n.end());
	NGT::ObjectID prev = 0;
	for (auto it = n.begin(); it != n.end();) {
	  if (prev == (*it).id) {
	    it = n.erase(it);
	    continue;
	  }
	  prev = (*it).id;
	  it++;
	}
	NGT::GraphNode tmp = n;
	n.swap(tmp);
      } catch (...) {
      }
    }
#endif
  }

  static void refineANNG(NGT::Index &index, bool unlog, float epsilon = 0.1, float accuracy = 0.0, int noOfEdges = 0, int exploreEdgeSize = INT_MIN, size_t batchSize = 10000) {
    NGT::StdOstreamRedirector redirector(unlog);
    redirector.begin();
    try {
      refineANNG(index, epsilon, accuracy, noOfEdges, exploreEdgeSize, batchSize);
    } catch (NGT::Exception &err) {
      redirector.end();
      throw(err);
    }
  }

  static void refineANNG(NGT::Index &index, float epsilon = 0.1, float accuracy = 0.0, int noOfEdges = 0, int exploreEdgeSize = INT_MIN, size_t batchSize = 10000) {
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
    NGTThrowException("GraphReconstructor::refineANNG: Not implemented for the shared memory option.");
#else
    auto prop = static_cast<GraphIndex&>(index.getIndex()).getGraphProperty();
    NGT::ObjectRepository &objectRepository = index.getObjectSpace().getRepository();
    NGT::GraphIndex &graphIndex = static_cast<GraphIndex&>(index.getIndex());
    size_t nOfObjects = objectRepository.size();
    bool error = false;
    std::string errorMessage;

    size_t noOfSearchedEdges = noOfEdges < 0 ? -noOfEdges : (noOfEdges > prop.edgeSizeForCreation ? noOfEdges : prop.edgeSizeForCreation);
    noOfSearchedEdges++;
    for (size_t bid = 1; bid < nOfObjects; bid += batchSize) {
      NGT::ObjectDistances results[batchSize];
      // search
#pragma omp parallel for
      for (size_t idx = 0; idx < batchSize; idx++) {
	size_t id = bid + idx;
	if (id % 100000 == 0) {
	  std::cerr << "# of processed objects=" << id << std::endl;
	}
	if (objectRepository.isEmpty(id)) {
	  continue;
	}
	NGT::SearchContainer searchContainer(*objectRepository.get(id));
	searchContainer.setResults(&results[idx]);
	assert(prop.edgeSizeForCreation > 0);
	searchContainer.setSize(noOfSearchedEdges);
	if (accuracy > 0.0) {
          searchContainer.setExpectedAccuracy(accuracy);
        } else {
	  searchContainer.setEpsilon(epsilon);
        }
	if (exploreEdgeSize != INT_MIN) {
          searchContainer.setEdgeSize(exploreEdgeSize);
        }
	if (!error) {
          try {
            index.search(searchContainer);
          } catch (NGT::Exception &err) {
#pragma omp critical
            {
      	      error = true;
	      errorMessage = err.what();
            }
          }
        }
      }
      if (error) {
        std::stringstream msg;
	msg << "GraphReconstructor::refineANNG: " << errorMessage;
        NGTThrowException(msg);
      }
      // outgoing edges
#pragma omp parallel for
      for (size_t idx = 0; idx < batchSize; idx++) {
	size_t id = bid + idx;
	if (objectRepository.isEmpty(id)) {
	  continue;
	}
	NGT::GraphNode &node = *graphIndex.getNode(id);
	for (auto i = results[idx].begin(); i != results[idx].end(); ++i) {
	  if ((*i).id != id) {
	    node.push_back(*i);
	  }
	}
	std::sort(node.begin(), node.end());
	// dedupe
	ObjectID prev = 0;
	for (GraphNode::iterator ni = node.begin(); ni != node.end();) {
	  if (prev == (*ni).id) {
	    ni = node.erase(ni);
	    continue;
	  }
	  prev = (*ni).id;
	  ni++;
	}
      }
      // incomming edges
      if (noOfEdges != 0) {
	continue;
      }
      for (size_t idx = 0; idx < batchSize; idx++) {
	size_t id = bid + idx;
	if (id % 10000 == 0) {
	  std::cerr << "# of processed objects=" << id << std::endl;
	}
	for (auto i = results[idx].begin(); i != results[idx].end(); ++i) {
	  if ((*i).id != id) {
	    NGT::GraphNode &node = *graphIndex.getNode((*i).id);
	    graphIndex.addEdge(node, id, (*i).distance, false);
	  }
	}
      }
    }
    if (noOfEdges > 0) {
      // prune to build knng
      size_t  nedges = noOfEdges < 0 ? -noOfEdges : noOfEdges;
#pragma omp parallel for
      for (ObjectID id = 1; id < nOfObjects; ++id) {
        if (objectRepository.isEmpty(id)) {
	  continue;
        }
	NGT::GraphNode &node = *graphIndex.getNode(id);
	if (node.size() > nedges) {
	  node.resize(nedges);
        }
      }
    }
#endif // defined(NGT_SHARED_MEMORY_ALLOCATOR)
  }
};

}; // NGT
