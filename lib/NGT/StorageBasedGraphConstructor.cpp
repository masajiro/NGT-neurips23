//
// Copyright (C) 2023 Yahoo Japan Corporation
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

#include "NGT/StorageBasedGraphConstructor.h"

void NGT::StorageBasedGraphConstructor::extractInsertionOrder(NGT::Index &index, size_t nOfThreads, size_t searchSize, float epsilon, bool indegreeOrder) {

  NGT::Timer timer;
  timer.start();
  std::cerr << "extractInsertionOrder" << std::endl;
  std::cerr << "VM size=" << NGT::Common::getProcessVmSizeStr() << std::endl;
  std::cerr << "Peak VM size=" << NGT::Common::getProcessVmPeakStr() << std::endl;
  std::cerr << "searching..." << std::endl;

  omp_set_num_threads(nOfThreads);

  std::cerr << "search size=" << searchSize << std::endl;
  if (index.getObjectRepositorySize() == 0) {
    std::stringstream msg;
    msg << "extractInsertionOrder: Error! No objects or no graph. " << index.getPath() << std::endl;
    NGTThrowException(msg);
  }
  std::vector<uint32_t> counter(nOfThreads);
  std::vector<uint32_t> indegrees[nOfThreads];
  for (size_t tidx = 0; tidx < nOfThreads; tidx++) {
    indegrees[tidx].resize(index.getObjectRepositorySize());
  }
  std::vector<std::pair<float, uint32_t>> length;
  length.resize(index.getObjectRepositorySize());
#pragma omp parallel for
  for (NGT::ObjectID query = 1; query < index.getObjectRepositorySize(); query++) {
    auto thdID = omp_get_thread_num();
    counter[thdID]++;
    if (query % 100000 == 0) {
      size_t n = 0;
      for (auto &c : counter) n += c;
      timer.stop();
      std::cerr << "# of the processed objects=" << n
		<< " VM size=" << NGT::Common::getProcessVmSizeStr()
		<< " Peak VM size=" << NGT::Common::getProcessVmPeakStr() 
		<< " Time=" << timer << std::endl;
      timer.restart();
    }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
    NGT::Object *object = index.getObjectSpace().allocateObject(*index.getObjectSpace().getRepository().get(query));
#else
    NGT::Object *object = index.getObjectSpace().getRepository().get(query);
#endif
    {
      NGT::SearchContainer sc(*object);
      NGT::ObjectDistances objects;
      sc.setResults(&objects);
      sc.setSize(searchSize);
      sc.setEpsilon(epsilon);
      sc.setEdgeSize(-2);
      NGT::Timer timer;
      try {
	timer.start();
	index.search(sc);
	timer.stop();
      } catch (NGT::Exception &err) {}
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      index.getObjectSpace().deleteObject(object);
#endif
      float len = 0.0;
      size_t count = 0;
      if (objects.size() == 0) {
	std::stringstream msg;
	msg << "extractInsertionOrder: Error! # of the searched objects is zero. " << index.getPath() << std::endl;
	NGTThrowException(msg);
      }
      for (size_t i = 0; i < objects.size(); i++) {
	if (query == objects[i].id) continue;
	len += objects[i].distance;
	count++;
	if (objects[i].id >= indegrees[thdID].size()) {
	  std::cerr << "too large. " << objects[i].id << ":" << indegrees[thdID].size() << std::endl;
	  exit(1);
	}
	indegrees[thdID][objects[i].id]++;
      }
      length[query].first = len / count;
      length[query].second = query;
    }
  }

  std::sort(length.begin(), length.end());

  size_t max = 0;
  for (NGT::ObjectID id = 1; id < index.getObjectRepositorySize(); id++) {
    for (size_t tidx = 1; tidx < nOfThreads; tidx++) {
      indegrees[0][id] += indegrees[tidx][id];
    }
    if (indegrees[0][id] > max) max = indegrees[0][id];
  }
  std::cerr << "max=" << max << std::endl;
  if (indegreeOrder) {
    std::vector<std::vector<uint32_t>> sortedIndegrees;
    sortedIndegrees.resize(max + 1);
    for (uint32_t oid = 1; oid < indegrees[0].size(); oid++) {

      sortedIndegrees[indegrees[0][oid]].push_back(oid);
    }
    {
      std::string path = index.getPath() + "/iord";
      size_t c = 0;
      std::ofstream of(path);
      if (!of) {
	std::stringstream msg;
	msg << "extractInsertionOrder: Cannot open the file. " << path << std::endl;
	NGTThrowException(msg);
      }
      for (uint32_t ind = 0; ind < sortedIndegrees.size(); ind++) {
	c += ind * sortedIndegrees[ind].size();
	if (sortedIndegrees[ind].size() != 0) {
	  for (auto &id: sortedIndegrees[ind]) {
	    of << id << "\t" << ind << std::endl;
	  }
	}
      }
      std::cerr << "total number of the incoming edges=" << c << ":" << (searchSize - 1) * (index.getObjectRepositorySize() - 1) << std::endl;
    }
  } else {
    std::string path = index.getPath() + "/iord";
    std::ofstream of(path);
    if (!of) {
      std::stringstream msg;
      msg << "extractInsertionOrder: Cannot open the file. " << path << std::endl;
      NGTThrowException(msg);
    }
    for (NGT::ObjectID id = index.getObjectRepositorySize() - 1; id != 0; id--) {
      of << length[id].second << "\t" << length[id].first << "\t" << indegrees[0][id] << std::endl;
    }
  }

}

void NGT::StorageBasedGraphConstructor::constructANNG(NGT::Index &index, size_t nOfThreads, size_t nOfGroups, size_t nOfMembers, size_t nOfSearchedObjects, std::vector<std::vector<std::string>> &paths, size_t outdegree, size_t indegree, bool verbose) {

  double totalTime	= 0;
  NGT::Timer timer;
  timer.start();

  paths.resize(nOfThreads);
  for (size_t t = 0; t < nOfThreads; t++) {
    paths[t].resize(nOfGroups);
    for (size_t d = 0; d < nOfGroups; d++) {
      std::stringstream path;
      path << index.getPath() << "/grp_" << t << "_" << d;
      paths[t][d] = path.str();
    }
  }

  std::string iordPath = index.getPath() + "/iord";
  std::ifstream iord(iordPath);
  std::vector<uint32_t> ids;
  auto validIds = false;
  if (iord) {
    std::cerr << "insertion order is valid." << std::endl;
    ids.reserve(index.getObjectRepositorySize());
    std::string line;
    while (true) {
      getline(iord, line);
      if (iord.eof()) break;
      std::vector<std::string> tokens;
      NGT::Common::tokenize(line, tokens, "\t, ");
      ids.push_back(NGT::Common::strtol(tokens[0]));
    }
    validIds = true;
  }

  auto &gtindex = static_cast<NGT::GraphAndTreeIndex &>(index.getIndex());
  auto batchSize = gtindex.GraphIndex::NeighborhoodGraph::property.batchSizeForCreation;
  std::cerr << "batchSize=" << batchSize << std::endl;
  size_t nOfEdges = gtindex.GraphIndex::NeighborhoodGraph::property.edgeSizeForCreation;
  std::cerr << "# of edges=" << nOfEdges << std::endl;
  auto nOfWrittenEdges = std::max(indegree, outdegree);

  std::ofstream ofs[nOfThreads][nOfGroups];
  for (size_t t = 0; t < nOfThreads; t++) {
    for (size_t d = 0; d < nOfGroups; d++) {
      ofs[t][d].open(paths[t][d], std::ios_base::out | std::ios::binary);
    }
  }
  std::cerr << "open is finished." << std::endl;
  std::cerr << "VM size=" << NGT::Common::getProcessVmSizeStr()
	    << " Peak VM size=" << NGT::Common::getProcessVmPeakStr() << std::endl;

  std::cerr << "extacting queries is finished." << std::endl;
  omp_set_num_threads(nOfThreads);
  std::cerr << "# of threads=" << nOfThreads << std::endl;
  std::cerr << "searching..." << std::endl;

  size_t nOfSeedNodes = 20;
  std::cerr << "insert the first " << nOfSeedNodes << " objects." << std::endl;
  for (NGT::ObjectID seedQuery = 1; seedQuery <= nOfSeedNodes; seedQuery++) {
    auto query = seedQuery;
    if (validIds) {
      query = ids[seedQuery - 1];
    }
    gtindex.insert(query);
  }
  std::vector<uint32_t> counter(nOfThreads);
  std::cerr << "insert the following objects." << std::endl;
  size_t notFoundCount = 0;
  for (NGT::ObjectID batchQuery = 1; batchQuery < index.getObjectRepositorySize(); batchQuery += batchSize) {
    NGT::ObjectDistances objects[batchSize];
#pragma omp parallel for
    for (int batch = 0; batch < batchSize; batch++) {
      auto query = batchQuery + batch;
      if (validIds) {
	query = ids[query - 1];
      }
      if (query >= index.getObjectRepositorySize()) {
	continue;
      }
      auto thdID = omp_get_thread_num();
      counter[thdID]++;
      if (query % 100000 == 0) {
	size_t n = 0;
	for (auto &c : counter) n += c;
	timer.stop();
	std::cerr << "# of the processed objects=" << n
		  << " VM size=" << NGT::Common::getProcessVmSizeStr()
		  << " Peak VM size=" << NGT::Common::getProcessVmPeakStr() 
		  << " Time=" << timer << std::endl;
	timer.restart();
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      NGT::Object *object = index.getObjectSpace().allocateObject(*index.getObjectSpace().getRepository().get(query));
#else
      NGT::Object *object = index.getObjectSpace().getRepository().get(query);
#endif
      {
	NGT::SearchContainer sc(*object);
	sc.setResults(&objects[batch]);
	sc.setSize(nOfSearchedObjects);
	sc.setRadius(std::numeric_limits<float>::max());
	sc.setEpsilon(gtindex.GraphIndex::NeighborhoodGraph::property.insertionRadiusCoefficient - 1.0);
	sc.setEdgeSize(-1);
	NGT::Timer searchTimer;
	try {
	  searchTimer.start();
	  index.search(sc);
	  searchTimer.stop();
	} catch (NGT::Exception &err) {}
	totalTime += searchTimer.time;
	if (verbose) {
	  std::cerr << "Query No." << query << std::endl;
	  std::cerr << "Rank\tID\tDistance" << std::endl;
	  for (size_t i = 0; i < objects[batch].size(); i++) {
	    std::cerr << i + 1 << "\t" << objects[batch][i].id << "\t";
	    std::cerr << objects[batch][i].distance << std::endl;
	  }
	  std::cerr << "Query Time= " << searchTimer.time << " (sec), " << searchTimer.time * 1000.0 << " (msec)" << std::endl;
	}
      }
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      index.getObjectSpace().deleteObject(object);
#endif
      uint32_t id = query;
      uint32_t gid = (query - 1) / nOfMembers;
      gid = gid >= nOfGroups ? nOfGroups - 1 : gid;
      uint32_t nid = id | 0x80000000;
      uint32_t count = 0;
      for (size_t i = 0; i < objects[batch].size() && count != nOfWrittenEdges; i++) {
	uint32_t dstid = objects[batch][i].id;
	if (dstid == id) continue;
	count++;
      }
      ofs[thdID][gid].write(reinterpret_cast<char*>(&nid), sizeof(nid));
      uint32_t len = count;
      ofs[thdID][gid].write(reinterpret_cast<char*>(&len), sizeof(len));
      count = 0;
      for (size_t i = 0; i < objects[batch].size() && count != len; i++) {
	uint32_t dstid = objects[batch][i].id;
	if (dstid == id) continue;
	float dstd = objects[batch][i].distance;
	ofs[thdID][gid].write(reinterpret_cast<char*>(&dstid), sizeof(dstid));
	ofs[thdID][gid].write(reinterpret_cast<char*>(&dstd), sizeof(dstd));
	count++;
      }
      bool found = false;
      count = 0;
      for (size_t i = 0; i < objects[batch].size() && count != nOfWrittenEdges; i++) {
	uint32_t dstid = objects[batch][i].id;
	if (dstid == id) {
	  found = true;
	  continue;
	}
	float dstd = objects[batch][i].distance;
	uint32_t dstgid = (dstid - 1) / nOfMembers;
	if (dstgid >= nOfGroups) {
	  dstgid = nOfGroups - 1;
	}
	ofs[thdID][dstgid].write(reinterpret_cast<char*>(&dstid), sizeof(dstid));
	ofs[thdID][dstgid].write(reinterpret_cast<char*>(&id), sizeof(id));
	ofs[thdID][dstgid].write(reinterpret_cast<char*>(&dstd), sizeof(dstd));
	count++;
      }
      if (!found) {
	notFoundCount++;
      }
    }
    for (int batch = 0; batch < batchSize; batch++) {
      auto query = batchQuery + batch;
      if (query <= nOfSeedNodes) {
	continue;
      }
      if (validIds) {
	query = ids[query - 1];
      }
      if (query >= index.getObjectRepositorySize()) {
	continue;
      }
      if (nOfEdges != 0 && objects[batch].size() > nOfEdges) {
	objects[batch].resize(nOfEdges);
      }
      gtindex.GraphIndex::NeighborhoodGraph::insertNode(query, objects[batch]);
#ifdef NGT_SHARED_MEMORY_ALLOCATOR
      NGT::Object *object = index.getObjectSpace().allocateObject(*index.getObjectSpace().getRepository().get(query));
      NGT::DVPTree::InsertContainer tiobj(*object, query);
      gtindex.DVPTree::insert(tiobj);
      index.getObjectSpace().deleteObject(object);
#else
      NGT::Object *object = index.getObjectSpace().getRepository().get(query);
      NGT::DVPTree::InsertContainer tiobj(*object, query);
      gtindex.DVPTree::insert(tiobj);
#endif
    }
  }
  if (verbose) {
    std::cerr << "Average Query Time= " << totalTime / (double)(index.getObjectRepositorySize() - 1)  << " (sec), "
	      << totalTime * 1000.0 / (double)(index.getObjectRepositorySize() - 1) << " (msec), ("
	      << totalTime << "/" << (index.getObjectRepositorySize() - 1) << ")" << std::endl;
  }
  for (size_t t = 0; t < nOfThreads; t++) {
    for (size_t d = 0; d < nOfGroups; d++) {
      uint32_t id = 0;
      ofs[t][d].write(reinterpret_cast<char*>(&id), sizeof(id));
    }
  }
  if (notFoundCount != 0) {
    std::cerr << "Warning: not found the target node in the result. " << notFoundCount << std::endl;
  }

  index.save();

}

void NGT::StorageBasedGraphConstructor::constructONNGFromANNG(std::string &indexPath, size_t &nOfThreads, size_t nOfGroups, size_t nOfMembers, std::vector<std::vector<std::string>> &paths, size_t outdegree, size_t indegree) {

  std::cerr << "constructONNGFromANNG..." << std::endl;

  std::string suffix = "_onng";
  if (nOfThreads == 0) {
    std::cerr << "nOfThread=0." << std::endl;
    exit(1);
  }
  if (nOfThreads != paths.size()) {
    std::cerr << "# of threads is inconsistency. " << nOfThreads << ":" << paths.size() << std::endl;
  }
  if (nOfGroups != paths[0].size()) {
    std::cerr << "# of groups is inconsistency. " << nOfGroups << ":" << paths[0].size() << std::endl;
  }
  std::ofstream ofs[nOfGroups];
  for (size_t d = 0; d < nOfGroups; d++) {
    ofs[d].open(paths[0][d] + suffix, std::ios_base::out | std::ios::binary);
  }

  for (size_t gid = 0; gid < nOfGroups; gid++) {
    std::vector<NGT::ObjectDistances> nodes;
    nodes.reserve(nOfMembers);
    size_t from = gid * nOfMembers + 1;
    size_t to = (gid + 1) * nOfMembers;
    std::cerr << "range=" << from << ":" << to << std::endl;
    for (size_t thdID = 0; thdID < nOfThreads; thdID++) {
      std::ifstream ifs(paths[thdID][gid], std::ios_base::in | std::ios::binary);
      do {
	uint32_t id;
	ifs.read(reinterpret_cast<char*>(&id), sizeof(id));
	if (id == 0) {
	  break;
	}
	uint32_t len = 1;
	if ((id & 0x80000000) != 0) {
	  id &= 0x7FFFFFFF;
	  ifs.read(reinterpret_cast<char*>(&len), sizeof(len));
	}
	if (id < from || id > to) {
	  std::stringstream msg;
	  msg << "Fatal error: out of range. " << id << ":" << from << "-" << to;
	  NGTThrowException(msg);
	}
	size_t idx = id - 1  - gid * nOfMembers;
	if (idx >= nodes.size()) {
	  nodes.resize(idx + 1);
	}
	for (size_t i = 0; i < len; i++) {
	  uint32_t dstid;
	  float dstd;
	  ifs.read(reinterpret_cast<char*>(&dstid), sizeof(dstid));
	  ifs.read(reinterpret_cast<char*>(&dstd), sizeof(dstd));
	  if (dstid > nOfGroups * nOfMembers + 1) {
	    std::cerr << "dstid is too large. " << dstid << std::endl;
	  }
	  nodes[idx].push_back(NGT::ObjectDistance(static_cast<size_t>(dstid), dstd));
	}
      } while (!ifs.eof());
    }
#pragma omp parallel for
    for (size_t ni = 0; ni < nodes.size(); ni++) {
      std::sort(nodes[ni].begin(), nodes[ni].end());
    }
    for (size_t ni = 0; ni < nodes.size(); ni++) {
      auto &node = nodes[ni];
      uint32_t id = gid * nOfMembers + ni + 1;
      uint32_t ogid = (id - 1) / nOfMembers;
      ogid = ogid >= nOfGroups ? nOfGroups - 1 : ogid;
      if (gid != ogid) {
	std::cerr << "something strange." << std::endl;
      }
      uint32_t nid = id | 0x80000000;
      uint32_t count = 0;
      uint32_t previd = 0;
      for (size_t i = 0; i < node.size(); i++) {
	uint32_t dstid = node[i].id;
	if (previd == node[i].id || dstid == id) {
	  node[i].id = 0;
	  previd = node[i].id;
	  continue;
	}
	previd = node[i].id;
	if (count >= outdegree) continue;
	count++;
      }
      if (count < outdegree) {
	std::cerr << "warning: too few outgoing edges. " << count << ":" << node.size() << std::endl;
      }
      ofs[gid].write(reinterpret_cast<char*>(&nid), sizeof(nid));
      uint32_t nOfEdges = count;
      ofs[gid].write(reinterpret_cast<char*>(&nOfEdges), sizeof(nOfEdges));
      count = 0;
      for (size_t i = 0; i < node.size() && count != nOfEdges; i++) {
	if (node[i].id == 0) continue;
	uint32_t dstid = node[i].id;
	float dstd = node[i].distance;
	ofs[gid].write(reinterpret_cast<char*>(&dstid), sizeof(dstid));
	ofs[gid].write(reinterpret_cast<char*>(&dstd), sizeof(dstd));
	count++;
      }
      count = 0;
      for (size_t i = 0; i < node.size() && count != indegree; i++) {
	if (node[i].id == 0) continue;
	uint32_t dstid = node[i].id;
	float dstd = node[i].distance;
	uint32_t dstgid = (dstid - 1) / nOfMembers;
	if (dstgid >= nOfGroups) {
	  dstgid = nOfGroups - 1;
	  std::cerr << "warning: a strange situation. " << dstid << std::endl;
	}
	ofs[dstgid].write(reinterpret_cast<char*>(&dstid), sizeof(dstid));
	ofs[dstgid].write(reinterpret_cast<char*>(&id), sizeof(id));
	ofs[dstgid].write(reinterpret_cast<char*>(&dstd), sizeof(dstd));
	count++;
      }
    }
  }

  for (size_t gid = 0; gid < nOfGroups; gid++) {
    uint32_t id = 0;
    ofs[gid].write(reinterpret_cast<char*>(&id), sizeof(id));
  }

  for (size_t gid = 0; gid < nOfGroups; gid++) {
    for (size_t thdID = 0; thdID < nOfThreads; thdID++) {
      unlink(paths[thdID][gid].c_str());
    }
  }
  {
    paths.resize(1);
    for (size_t gid = 0; gid < nOfGroups; gid++) {
      paths[0][gid] = paths[0][gid] + suffix;
    }
  }
}

void NGT::StorageBasedGraphConstructor::aggregate(std::string &indexPath, size_t nOfThreads, size_t nOfGroups, size_t nOfMembers, std::vector<std::vector<std::string>> &paths) {

  std::cerr << "aggregating..." << std::endl;

  if (paths.size() == 0 || paths[0].size() == 0) {
    std::cerr << "paths are invalid. " << std::endl;
    exit(1);
  }
  for (size_t gid = 0; gid < nOfGroups; gid++) {
    std::vector<NGT::ObjectDistances> nodes;
    nodes.reserve(nOfMembers + nOfGroups);
    size_t from = gid * nOfMembers + 1;
    size_t to = (gid + 1) * nOfMembers;
    std::cerr << "range=" << from << ":" << to << std::endl;
    for (size_t thdID = 0; thdID < paths.size(); thdID++) {
      std::ifstream ifs(paths[thdID][gid], std::ios_base::in | std::ios::binary);
      do {
	uint32_t id;
	ifs.read(reinterpret_cast<char*>(&id), sizeof(id));
	if (id == 0) {
	  break;
	}
	uint32_t len = 1;
	if ((id & 0x80000000) != 0) {
	  id &= 0x7FFFFFFF;
	  ifs.read(reinterpret_cast<char*>(&len), sizeof(len));
	}
	if (id < from || id > to) {
	  std::cerr << "aggregate: Warning: out of range. " << id << ":" << from << "-" << to << std::endl;
	  float v;
	  ifs.read(reinterpret_cast<char*>(&v), sizeof(v));
	  continue;
	}
	size_t idx = id - 1  - gid * nOfMembers;
	if (idx >= nodes.size()) {
	  nodes.resize(idx + 1);
	}
	for (size_t i = 0; i < len; i++) {
	  uint32_t dstid;
	  float dstd;
	  ifs.read(reinterpret_cast<char*>(&dstid), sizeof(dstid));
	  ifs.read(reinterpret_cast<char*>(&dstd), sizeof(dstd));
	  if (dstid > nOfGroups * nOfMembers + 1) {
	    std::cerr << "dst node ID is too large. " << dstid << ":" << nOfGroups * nOfMembers << std::endl;
	    std::cerr << "     i=" << i << " ID=" << id << " gid=" << gid << " thdID=" << thdID << " len=" << len << std::endl;
	    std::cerr << "     " << ((dstid & 0xFFFF0000) >> 16) << ":" << (dstid & 0xFFFF) << std::endl;
	    exit(1);
	  }
	  nodes[idx].push_back(NGT::ObjectDistance(static_cast<size_t>(dstid), dstd));
	}
      } while (!ifs.eof());
    }
#pragma omp parallel for
    for (size_t ni = 0; ni < nodes.size(); ni++) {
      std::sort(nodes[ni].begin(), nodes[ni].end());
    }

    std::stringstream path;
    path << indexPath << "/grp_" << gid;
    std::ofstream of(path.str(), std::ios_base::out | std::ios::binary);
    for (size_t ni = 0; ni < nodes.size(); ni++) {
      auto &node = nodes[ni];
      uint32_t id = gid * nOfMembers + ni + 1;
      uint32_t count = 0;
      uint32_t previd = 0;
      for (size_t i = 0; i < node.size(); i++) {
	if (previd == node[i].id) {
	  node[i].id = 0;
	} else {
	  count++;
	}
	previd = node[i].id;
      }

      of.write(reinterpret_cast<char*>(&id), sizeof(id));
      if (node.size() == 0) {
	std::cerr << "warning: empty node! " << ni << ":" << id << std::endl;
      }
      uint32_t size = count;
      of.write(reinterpret_cast<char*>(&size), sizeof(size));
      float prevDistance = 0.0;
      for (auto edge : node) {
	uint32_t eid = edge.id;
	if (eid == 0) continue;
	float edistance = edge.distance;
	if (prevDistance > edistance) {
	  std::cerr << "Fatal: the distance order is invalid. " << id << ":" << eid << ":" << edistance << ":" << prevDistance << std::endl;
	}
	of.write(reinterpret_cast<char*>(&eid), sizeof(eid));
	of.write(reinterpret_cast<char*>(&edistance), sizeof(edistance));
      }
    }
  }
  for (size_t thdID = 0; thdID < paths.size(); thdID++) {
    for (size_t gid = 0; gid < paths[thdID].size(); gid++) {
      unlink(paths[thdID][gid].c_str());
    }
  }
}

void NGT::StorageBasedGraphConstructor::concatinate(std::string &indexPath, size_t nOfObjects, size_t nOfGroups, size_t nOfFinalEdges) {

  std::cerr << "concatinating..." << std::endl;

  std::stringstream path;
  path << indexPath << "/grp_";
  std::ofstream ofs(path.str(), std::ios_base::out | std::ios::binary);
  if (!ofs) {
    std::cerr << "cannot open. " << path.str() << std::endl;
  }
  uint64_t repoSize = nOfObjects + 1;
  ofs.write(reinterpret_cast<char*>(&repoSize), sizeof(repoSize));
  char flag = '-';
  ofs.write(reinterpret_cast<char*>(&flag), sizeof(flag));

  std::cerr << "# of groups=" << nOfGroups << std::endl;
  for (size_t gid = 0; gid < nOfGroups; gid++) {
    std::stringstream path;
    path << indexPath << "/grp_" << gid;
    std::ifstream ifs(path.str(), std::ios_base::in | std::ios::binary);
    std::cerr << "path=" << path.str() << std::endl;
    if (!ifs) {
      std::cerr << "cannot open. " << path.str() << std::endl;
      exit(1);
    }
    uint32_t prevID = 0;
    for (;;) {
      uint32_t id = 0;
      ifs.read(reinterpret_cast<char*>(&id), sizeof(id));
      if (prevID != 0 && prevID + 1 != id) {
	std::cerr << "invalid ID. " << prevID << ":" << id << std::endl;
	exit(1);
      }
      if (ifs.eof()) {
	break;
      }
      flag = '+';
      ofs.write(reinterpret_cast<char*>(&flag), sizeof(flag));
      uint32_t size;
      ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
      uint32_t esize = (nOfFinalEdges != 0 && size >= nOfFinalEdges) ? nOfFinalEdges : size;
      ofs.write(reinterpret_cast<char*>(&esize), sizeof(esize));
      for (size_t ei = 0; ei < size; ei++) {
	uint32_t eid;
	float edistance;
	ifs.read(reinterpret_cast<char*>(&eid), sizeof(eid));
	ifs.read(reinterpret_cast<char*>(&edistance), sizeof(edistance));
	if (eid >= repoSize) {
	  std::cerr << "Fatal error! edge ID is too large. " << eid << ":" << repoSize << std::endl;
	  exit(1);
	}
	if (nOfFinalEdges != 0 && ei >= nOfFinalEdges) continue;
	ofs.write(reinterpret_cast<char*>(&eid), sizeof(eid));
	ofs.write(reinterpret_cast<char*>(&edistance), sizeof(edistance));
      }
    }
    unlink(path.str().c_str());
  }

}

void NGT::StorageBasedGraphConstructor::construct(std::string &indexPath, std::string &baseIndexPath, char graphType, size_t nOfSearchedObjects, size_t outdegree, size_t indegree, size_t nOfFinalEdges, size_t searchSize, float epsilon, size_t nOfThreads, bool verbose) {

  StdOstreamRedirector redirector(!verbose);
  redirector.begin();

  try {
    if (!baseIndexPath.empty()) {
      std::string objPath = indexPath + "/obj";
      std::string iordPath = indexPath + "/iord";
      std::string baseObjPath;
      std::string baseIordPath;
      unlink(objPath.c_str());
      unlink(iordPath.c_str());
      if (baseIndexPath[0] != '/') {
	int i;
	for (i = 0; static_cast<size_t>(i) < baseIndexPath.size() && 
	       static_cast<size_t>(i) < indexPath.size(); i++) {
	  if (baseIndexPath[i] != indexPath[i]) break;
	}
	for (; i > 0; i--) {
	  if (baseIndexPath[i] == '/') {
	    i++;
	    break;
	  }
	}
	int cnt = 0;
	for (int j = indexPath.size() - 1; j >= i; j--) {
	  if (indexPath[j] == '/') cnt++;
	}
	cnt++;
	std::string p;
	for (int c = 0; c < cnt; c++) {
	  p += "../";
	}
	baseObjPath = p + baseIndexPath.substr(i) + "/obj";
	baseIordPath = p + baseIndexPath.substr(i) + "/iord";
      }
      symlink(baseObjPath.c_str(), objPath.c_str());
      symlink(baseIordPath.c_str(), iordPath.c_str());
    }

    if (nOfThreads == 0 || nOfThreads > static_cast<size_t>(omp_get_max_threads())) {
      nOfThreads = omp_get_max_threads();
    }

    size_t nOfGroups = 10;
    size_t nOfObjects;
    size_t nOfMembers;
    std::vector<std::vector<std::string>> paths;
    {
      NGT::Index index(indexPath, false);
      indexPath = index.getPath();
      nOfObjects = index.getObjectRepositorySize() - 1;
      nOfMembers = nOfObjects / nOfGroups;
      if (nOfObjects % nOfGroups != 0) {
	nOfMembers++;
      }
      if (verbose) {
	std::cerr << "index is open " << std::endl;
	std::cerr << "VM size=" << NGT::Common::getProcessVmSizeStr() << std::endl;
	std::cerr << "Peak VM size=" << NGT::Common::getProcessVmPeakStr() << std::endl;
	std::cerr << "# of searched objects=" << nOfSearchedObjects << std::endl;
	std::cerr << "outdegree=" << outdegree << std::endl;
	std::cerr << "indegree=" << indegree << std::endl;
	std::cerr << "# of objects=" << nOfObjects << std::endl;;
	std::cerr << "# of threads=" << nOfThreads << std::endl;
	std::cerr << "# of groups=" << nOfGroups << std::endl;
	std::cerr << "# of members=" << nOfMembers << std::endl;
      }

      if (graphType == '-') {
	NGT::StorageBasedGraphConstructor::extractInsertionOrder(index, nOfThreads, searchSize, epsilon);
	redirector.end();
	return;
      }

      NGT::StorageBasedGraphConstructor::constructANNG(index, nOfThreads, nOfGroups, nOfMembers, nOfSearchedObjects, paths, outdegree, indegree);
      if (verbose) {
	std::cerr << "search is finished." << std::endl;
	std::cerr << "VM size=" << NGT::Common::getProcessVmSizeStr() << std::endl;
	std::cerr << "Peak VM size=" << NGT::Common::getProcessVmPeakStr() << std::endl;
      }
    }
    if (graphType == 'o') {
      NGT::StorageBasedGraphConstructor::constructONNGFromANNG(indexPath, nOfThreads, nOfGroups, nOfMembers, paths, outdegree, indegree);
    } if (graphType == 'a') {
      NGT::StorageBasedGraphConstructor::constructONNGFromANNG(indexPath, nOfThreads, nOfGroups, nOfMembers, paths, 
							       std::numeric_limits<size_t>::max(), 
							       std::numeric_limits<size_t>::max());
    } else {
      std::cerr << "invalid graph type. " << graphType << std::endl;
    }
    NGT::StorageBasedGraphConstructor::aggregate(indexPath, nOfThreads, nOfGroups, nOfMembers, paths);
    if (verbose) {
      std::cerr << "aggregation is finished." << std::endl;
      std::cerr << "VM size=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "Peak VM size=" << NGT::Common::getProcessVmPeakStr() << std::endl;
    }
    NGT::StorageBasedGraphConstructor::concatinate(indexPath, nOfObjects, nOfGroups, nOfFinalEdges);
    if (verbose) {
      std::cerr << "concatination is finished." << std::endl;
      std::cerr << "VM size=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "Peak VM size=" << NGT::Common::getProcessVmPeakStr() << std::endl;
    }
    {
      std::string src = indexPath + "/grp_";
      std::string dst = indexPath + "/grp";
      rename(src.c_str(), dst.c_str());
    }
  } catch (NGT::Exception &err) {
    redirector.end();
    throw err;
  }
  redirector.end();
}
