//
// Created by Gabriele Canesi on 01/07/23.
//


#include "util.h"
#include <string>
#include <vector>


std::vector<std::string> util::buildAllSpacedPatterns(size_t k) {
    std::vector<std::string> result(k, std::string(k, '1'));
    for (size_t i = 0; i < k; ++i) {
        result[i][i] = '0';
    }

    return result;
}


std::vector<std::string> util::buildFirstSpacedPattern(size_t k) {
    std::vector<std::string> result(1, std::string(k, '1'));
    result[0][0] = '0';

    return result;
}

std::vector<std::string> util::buildRemainingSpacedPatterns(size_t k) {
    std::vector<std::string> result(k - 1, std::string(k, '1'));
    for (size_t i = 1; i < k; ++i) {
        result[i - 1][i] = '0';
    }

    return result;
}

void util::addToHashTable(uint64_t hash, size_t position, std::unordered_map<uint64_t, std::vector<size_t>> &hashTable) {
    auto iterator = hashTable.find(hash);
    if (iterator == hashTable.end()) {
        hashTable[hash] = std::vector<size_t>(1, position);
    } else {
        hashTable[hash].push_back(position);
    }
}