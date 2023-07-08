//
// Created by Gabriele Canesi on 07/07/23.
//

#include "MinHashInfo.h"
#include <ankerl/unordered_dense.h>
#include <nthash/nthash.hpp>
#include <random>
#include <vector>

MinHashInfo::MinHashInfo(const std::string &reference, const std::string &r, size_t k) {

    auto intervals = reference.length() / r.length();

    maps = std::vector<ankerl::unordered_dense::set<uint64_t, hash>>(intervals + 1, ankerl::unordered_dense::set<uint64_t, hash>());
    nthash::NtHash nth(r, 1, k);
    while (nth.roll()) {
        maps[0].insert(nth.hashes()[0]);
    }

    std::vector<std::vector<uint64_t>> seqKmers = {};
    seqKmers.reserve(intervals);

    size_t seq = 0;
    size_t index = 0;

    nthash::NtHash nth2(reference, 1, k);

    while (nth2.roll() && seq < intervals) {
        auto &present = maps[seq + 1];
        if (index <= r.length() - k) {
            present.insert(nth2.hashes()[0]);
        }
        ++index;
        if (r.length() - k + 1 == index) {
            ++seq;
            index = 0;
        }
    }


    int i = 0;for (auto &map : maps) {
        for (auto el : map) {
            allKMers.insert(el);
        }
        ++i;
    }

    std::cout << "Ended inserting" << std::endl;

}


void MinHashInfo::computeSignature() {
    signature = std::vector<std::vector<uint64_t>>(NUM_HASHES, std::vector<uint64_t>(maps.size(), 0));
    std::random_device rd;
    std::mt19937  g(rd());
    std::vector<size_t> values(allKMers.begin(), allKMers.end());


    for (size_t i = 0; i < NUM_HASHES; ++i) {
        std::shuffle(values.begin(), values.end(), g);
        size_t currentMap = 0;
        for (auto &map : maps) {
                auto mapIterator = map.find(values[i]);
                size_t j = 0;
                while (mapIterator == map.end()){
                    ++j;
                    mapIterator = map.find(values[(i + j) % values.size()]);
                }
                signature[i][currentMap] = j;
            ++currentMap;
        }
    }
}


std::vector<double> MinHashInfo::computeJaccard() {
    std::vector<double> result(signature[0].size() - 1, 0);
    for (int i = 1; i < signature[0].size(); ++i) {
        size_t matches = 0;
        for (auto hash = 0; hash < NUM_HASHES; ++hash) {
            if (signature[hash][0] == signature[hash][i]) {
                ++matches;
            }
        }
        result[i - 1] = (double) matches / (double) NUM_HASHES;
    }

    return result;
}