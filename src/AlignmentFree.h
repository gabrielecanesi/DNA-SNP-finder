//
// Created by Gabriele Canesi on 07/07/23.
//

#ifndef SPACEDSEEDS_ALIGNMENTFREE_H
#define SPACEDSEEDS_ALIGNMENTFREE_H

#include <vector>
#include <ankerl/unordered_dense.h>
#include <cstddef>
#include <unordered_set>
#include <set>
#include "BloomFilter.h"

class AlignmentFree {
    BloomFilter substringBloomFilter;
    std::vector<double> jaccards;
    std::vector<uint64_t> M_rKmers;

public:
    AlignmentFree(const std::string &reference, const std::string &r, size_t k, size_t size);
    const std::vector<uint64_t>& rKmers();
    const std::vector<double>& getJaccard();
    const BloomFilter& getSubstringBloomFilter();
};

#endif //SPACEDSEEDS_ALIGNMENTFREE_H
