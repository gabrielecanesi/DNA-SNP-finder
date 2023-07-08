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

constexpr uint64_t NUM_HASHES = 200;

struct hash {
    auto operator()(uint64_t const& x) const noexcept -> uint64_t {
        return x;
    }
};

class MinHashInfo {
    ankerl::unordered_dense::set<uint64_t, hash> allKMers;
    std::vector<ankerl::unordered_dense::set<uint64_t, hash>> maps;
    std::vector<std::vector<uint64_t>> signature;

public:
    MinHashInfo(const std::string &reference, const std::string &r, size_t k);
    void computeSignature();
    std::vector<double> computeJaccard();
};

#endif //SPACEDSEEDS_MINHASHINFO_H
