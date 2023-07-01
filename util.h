//
// Created by Gabriele Canesi on 01/07/23.
//

#ifndef SPACED_SEEDS_UTIL_H
#define SPACED_SEEDS_UTIL_H

#include <string>
#include <vector>

namespace util {
    std::vector<std::string> buildAllSpacedPatterns(size_t k);
    std::vector<std::string> buildFirstSpacedPattern(size_t k);
    std::vector<std::string> buildRemainingSpacedPatterns(size_t k);
    void addToHashTable(uint64_t hash, size_t position, std::unordered_map<uint64_t, std::vector<size_t>> &hashTable);
}


#endif //SPACED_SEEDS_UTIL_H
