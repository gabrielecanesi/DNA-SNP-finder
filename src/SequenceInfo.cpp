//
// Created by Gabriele Canesi on 01/07/23.
//

#include "SequenceInfo.h"
#include "util.h"
#include "nthash/nthash.hpp"

SequenceInfo::SequenceInfo(const char *sequence, const std::vector<std::string> &seeds, size_t length, bool buildHashTable) :
        M_sequence(sequence), firstSeedsHashes(), hashTable(), M_k(seeds[0].length()), sequence_length(length){

    nthash::SeedNtHash nth(sequence, sequence_length, seeds, 1, M_k);
    int j = 1;
    while (nth.roll()) {
        if (buildHashTable){
            for (size_t i = 0; i < k(); ++i) {
               util::addToHashTable(nth.hashes()[i], nth.get_pos(), hashTable);
            }

            ++j;
        }
        firstSeedsHashes.push_back(nth.hashes()[0]);
    }
}

SequenceInfo SequenceInfo::buildForReference(const char *sequence, size_t k, size_t length) {
    auto pattern = util::buildAllSpacedPatterns(k);
    return {sequence, pattern, length};
}

SequenceInfo SequenceInfo::buildForSubstring(const char *sequence, size_t k, size_t length) {
    auto pattern = util::buildFirstSpacedPattern(k);
    SequenceInfo info = {sequence, pattern, length, false};
    info.extendLastSeeds();
    return info;
}

void SequenceInfo::extendLastSeeds() {
    auto pattern = util::buildRemainingSpacedPatterns(k());
    nthash::SeedNtHash nth(M_sequence, pattern, 1, k(), sequence_length - k());
    nth.roll();
    for (int i = 0; i < k() - 1; ++i) {
        firstSeedsHashes.push_back(nth.hashes()[i]);
    }
}

const char * SequenceInfo::sequence() {
    return M_sequence;
}

const std::list<size_t>& SequenceInfo::positionsForHash(uint64_t hash) const {
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