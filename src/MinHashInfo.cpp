//
// Created by Gabriele Canesi on 07/07/23.
//

#include "MinHashInfo.h"
#include <ankerl/unordered_dense.h>
#include <nthash/nthash.hpp>
#include <random>
#include <vector>

MinHashInfo::MinHashInfo(const std::string &reference, const std::string &r, size_t k) : substringBloomFilter(10000000) {


    auto intervals = reference.length() / r.length();
    jaccards = std::vector<double>(intervals, 0.0);

    nthash::NtHash nth(r, 1, k);
    while (nth.roll()) {
        substringBloomFilter.insert(nth.hashes()[0]);
    }

    nthash::NtHash nth2(reference, 1, k);

    for (size_t seq = 0; seq < intervals; ++seq) {
        size_t matches = 0;
        for (size_t j = 0; j <= r.length() - k; ++j) {
            nth2.roll();
            if (substringBloomFilter.contains(nth2.hashes()[0])) {
                ++matches;
            }
        }
        jaccards[seq] = (double) matches / (double) (r.length() - k);

        for (size_t j = 0; j < k - 1; ++j) {
            nth2.roll();
        }
    }

    std::cout << "Ended insertion and computation of Jaccard metric" << std::endl;
}


void MinHashInfo::computeSignature() {

}


const std::vector<double>& MinHashInfo::computeJaccard() {
    return jaccards;
}



const std::vector<size_t>& MinHashInfo::candidates() {
    return M_candidates;
}