//
// Created by Gabriele Canesi on 01/07/23.
//

#include "SequenceInfo.h"
#include "util.h"
#include "nthash/nthash.hpp"
#include "BloomFilter.h"
#include "util.h"

SequenceInfo::SequenceInfo(const char *sequence, const std::vector<std::string> &seeds, size_t length,  const ankerl::unordered_dense::set<uint64_t>& toCompare, bool buildHashTable) :
        M_sequence(sequence), firstSeedsHashes(), hashTable(), M_k(seeds[0].length()), sequence_length(length), filter(0){

    ankerl::unordered_dense::map<uint64_t, size_t> frequencies;
    nthash::SeedNtHash nth(sequence, sequence_length, seeds, 1, M_k);

    while (nth.roll()) {
        for (size_t i = 0; i < k(); ++i) {
            if (toCompare.contains(nth.hashes()[i])) {
                util::addToHashTable(nth.hashes()[i], nth.get_pos(), hashTable);

                /*auto iterator = frequencies.find(nth.hashes()[i]);
                if (iterator == frequencies.end()) {
                    frequencies.insert({nth.hashes()[i], 1});
                } else {
                    ++iterator->second;
                }*/
           }
        }
        firstSeedsHashes.push_back(nth.hashes()[0]);
    }

    /*for (auto pair : frequencies) {
        orderedRHashes.push_back(pair);
    }

    std::sort(orderedRHashes.begin(), orderedRHashes.end(), [](const auto &x, const auto &y){
        return x.second < y.second;
    });*/
}

SequenceInfo::SequenceInfo(const char *sequence, const std::vector<std::string> &seeds, size_t length, bool buildHashTable) :
        M_sequence(sequence), firstSeedsHashes(), hashTable(), M_k(seeds[0].length()), sequence_length(length), filter(){

    nthash::SeedNtHash nth(sequence, sequence_length, seeds, 1, M_k, 0);
    for (size_t i = 0; i <= length - k(); ++i) {
        size_t limit = i < length - k() ? 1 : k();
        nth.roll();
        for (int j = 0; j < limit; ++j) {
            filter.insert(nth.hashes()[j]);
            firstSeedsHashes.push_back(nth.hashes()[j]);
            util::addToHashTable(nth.hashes()[j], nth.get_pos(), hashTable);
        }
    }
}

SequenceInfo SequenceInfo::buildForReference(const char *sequence, size_t k, size_t length, const ankerl::unordered_dense::set<uint64_t>& rFilter) {
    auto pattern = util::buildAllSpacedPatterns(k);
    return {sequence, pattern, length, rFilter, true};
}

SequenceInfo SequenceInfo::buildForSubstring(const char *sequence, size_t k, size_t length) {
    auto pattern = util::buildAllSpacedPatterns(k);
    return {sequence, pattern, length, false};
}


const char * SequenceInfo::sequence() {
    return M_sequence;
}

const std::vector<size_t>& SequenceInfo::positionsForHash(uint64_t hash) const {
    auto iterator = hashTable.find(hash);
    if (iterator == hashTable.end()) {
        return nullVector;
    }

    return iterator->second;
}

uint64_t SequenceInfo::hashAtPosition(size_t position) const {
    return firstSeedsHashes[position];
}

size_t SequenceInfo::k() const {
    return M_k;
}

size_t SequenceInfo::sequenceLength() const {
    return sequence_length;
}

const ankerl::unordered_dense::set<uint64_t>& SequenceInfo::getFilter() const {
    return filter;
}

const std::vector<std::pair<uint64_t, size_t>>& SequenceInfo::getOrderedRHashes() const {
    return orderedRHashes;
}