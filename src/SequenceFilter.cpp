//
// Created by Gabriele Canesi on 07/07/23.
//

#include "SequenceFilter.h"
#include <ankerl/unordered_dense.h>
#include <nthash/nthash.hpp>
#include <vector>
#include <cmath>

SequenceFilter::SequenceFilter(const std::string &reference, const std::string &r, size_t k, size_t size, double bloomFilterThreshold)  {


    bloom_parameters parameters;
    parameters.projected_element_count = std::min(r.length(), (size_t) std::pow(4, k));
    parameters.false_positive_probability = bloomFilterThreshold;
    parameters.compute_optimal_parameters();
    substringBloomFilter = bloom_filter(parameters);

    size_t intervals = reference.length() / size;
    similarities = std::vector<double>(intervals, 0.0);

    nthash::NtHash nth(r, 1, k);
    M_rKmers = std::vector<uint64_t>(r.length() - k + 1, 0);
    size_t i = 0;
    while (nth.roll()) {
        substringBloomFilter.insert(nth.hashes()[0]);
        M_rKmers[i] = nth.hashes()[0];
        ++i;
    }

    nthash::NtHash nth2(reference, 1, k);

    for (size_t seq = 0; seq < intervals; ++seq) {
        size_t matches = 0;
        for (size_t j = 0; j <= size - k; ++j) {
            nth2.roll();
            if (substringBloomFilter.contains(nth2.hashes()[0])) {
                ++matches;
            }
        }
        similarities[seq] = (double) matches / (double) (size);

        for (size_t j = 0; j < k - 1; ++j) {
            nth2.roll();
        }
    }
}


const std::vector<double>& SequenceFilter::getSimilarities() {
    return similarities;
}

const std::vector<uint64_t>& SequenceFilter::rKmers() {
    return M_rKmers;
}


const bloom_filter& SequenceFilter::getSubstringBloomFilter() {
    return substringBloomFilter;
}
