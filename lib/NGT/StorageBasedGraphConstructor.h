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

#pragma once

#include	<unordered_map>
#include	<unordered_set>
#include	<list>

#ifdef _OPENMP
#include	<omp.h>
#else
#warning "*** OMP is *NOT* available! ***"
#endif

#include "NGT/Index.h"

namespace NGT {
  class StorageBasedGraphConstructor {
  public:
    StorageBasedGraphConstructor(){}
    static void extractInsertionOrder(NGT::Index &index, size_t nOfThreads, size_t searchSize = 100, float epsilon = 0.01, bool indegreeOrder = false);
    static void constructANNG(NGT::Index &index, size_t nOfThreads, size_t nOfGroups, size_t nOfMembers, size_t nOfSearchedObjects, std::vector<std::vector<std::string>> &paths, size_t outdegree, size_t indegree, bool verbose = false);
    static void constructONNGFromANNG(std::string &indexPath, size_t &nOfThreads, size_t nOfGroups, size_t nOfMembers, std::vector<std::vector<std::string>> &paths, size_t outdegree, size_t indegree);
    static void aggregate(std::string &indexPath, size_t nOfThreads, size_t nOfGroups, size_t nOfMembers, std::vector<std::vector<std::string>> &paths);
    static void concatinate(std::string &indexPath, size_t nOfObjects, size_t nOfGroups, size_t nOfFinalEdges);
    static void construct(std::string &indexPath, std::string &baseIndexPath, char graphType, size_t nOfSearchedObjects, size_t outdegree, size_t indegree, size_t nOfFinalEdges, size_t searchSize, float epsilon, size_t nOfThreads = 0, bool verbose = false);
  };

}
