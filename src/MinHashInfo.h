//
// Created by Gabriele Canesi on 07/07/23.
//

#ifndef SPACEDSEEDS_MINHASHINFO_H
#define SPACEDSEEDS_MINHASHINFO_H

#include <vector>
#include <ankerl/unordered_dense.h>
#include <cstddef>
#include <unordered_set>
#include <set>
#include "BloomFilter.h"

constexpr uint64_t NUM_HASHES = 100;

struct hash {
    auto operator()(uint64_t const& x) const noexcept -> uint64_t {
        return x;
    }
};

class MinHashInfo {
    ankerl::unordered_dense::set<uint64_t, hash> allKMers;
    std::vector<ankerl::unordered_dense::set<uint64_t, hash>> maps;
    std::vector<std::vector<uint64_t>> signature;
    std::vector<size_t> M_candidates;
    BloomFilter substringBloomFilter;
    std::vector<double> jaccards;

public:
    MinHashInfo(const std::string &reference, const std::string &r, size_t k);
    void computeSignature();
    const std::vector<double>& computeJaccard();
    const std::vector<size_t>& candidates();
};

#endif //SPACEDSEEDS_MINHASHINFO_H
