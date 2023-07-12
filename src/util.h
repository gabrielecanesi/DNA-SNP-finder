//
// Created by Gabriele Canesi on 01/07/23.
//

#ifndef SPACED_SEEDS_UTIL_H
#define SPACED_SEEDS_UTIL_H

#include <string>
#include <vector>
#include <list>
#include <ankerl/unordered_dense.h>
#include <iostream>

template <typename  K, typename V>
std::ostream& operator<<(std::ostream &out, const std::pair<K, V> &pair) {
    out << "(" << pair.first << ", " << pair.second << ")";
    return out;
}

template <typename  T>
std::ostream& operator<<(std::ostream &out, const std::vector<T> &vector) {
    out << "[";
    for (size_t i = 0; i < vector.size() && vector.size() > 0; ++i) {
        out << vector[i] << ", ";
    }
    if (vector.size() > 0) {
        out << vector[vector.size() - 1];
    }
    out << "]";

    return out;
}



namespace util {
    std::vector<std::string> buildAllSpacedPatterns(size_t k);
    std::vector<std::string> buildFirstSpacedPattern(size_t k);
    std::vector<std::string> buildRemainingSpacedPatterns(size_t k);
    inline void addToHashTable(uint64_t hash, size_t position, ankerl::unordered_dense::map<uint64_t, std::vector<size_t>> &hashTable) {

        auto iterator = hashTable.find(hash);
        if (iterator == hashTable.end()) {
            hashTable.insert({hash, {position}});
        } else {
            iterator->second.push_back(position);
        }
    }
    std::shared_ptr<std::string> readFromFASTA(const std::string &path, bool skipNonACGT);
}


#endif //SPACED_SEEDS_UTIL_H
