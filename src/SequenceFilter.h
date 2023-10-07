//
// Created by Gabriele Canesi on 07/07/23.
//

#ifndef SPACEDSEEDS_SEQUENCEFILTER_H
#define SPACEDSEEDS_SEQUENCEFILTER_H

#include <vector>
#include <ankerl/unordered_dense.h>
#include <cstddef>
#include <unordered_set>
#include <set>
#include <bloom/bloom_filter.hpp>

class SequenceFilter {
    bloom_filter substringBloomFilter;
    std::vector<double> similarities;
    std::vector<uint64_t> M_rKmers;

public:
    SequenceFilter(const std::string &reference, const std::string &r, size_t k, size_t size,
                   double bloomFilterThreshold);
    const std::vector<uint64_t>& rKmers();
    const std::vector<double>& getSimilarities();
    const bloom_filter& getSubstringBloomFilter();
};

#endif //SPACEDSEEDS_SEQUENCEFILTER_H
