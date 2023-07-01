//
// Created by Gabriele Canesi on 01/07/23.
//

#include "SequenceInfo.h"
#include "util.h"
#include <nthash/nthash.hpp>

SequenceInfo::SequenceInfo(const std::shared_ptr<std::string> &sequence, const std::vector<std::string> &seeds, bool buildHashTable) :
        M_sequence(sequence), firstSeedsHashes(), hashTable(), M_k(seeds[0].length()){

    nthash::SeedNtHash nth(*sequence, seeds, 1, M_k);
    while (nth.roll()) {
        if (buildHashTable){
            for (size_t i = 0; i < k(); ++i) {
                util::addToHashTable(nth.hashes()[i], nth.get_pos(), hashTable);
            }
        }
        firstSeedsHashes.push_back(nth.hashes()[0]);
    }
}

SequenceInfo SequenceInfo::buildForReference(const std::shared_ptr<std::string> &sequence, size_t k) {
    auto pattern = util::buildAllSpacedPatterns(k);
    return {sequence, pattern};
}

SequenceInfo SequenceInfo::buildForSubstring(const std::shared_ptr<std::string> &sequence, size_t k) {
    auto pattern = util::buildFirstSpacedPattern(k);
    SequenceInfo info = {sequence, pattern, false};
    info.extendLastSeeds();
    return info;
}

void SequenceInfo::extendLastSeeds() {
    auto pattern = util::buildRemainingSpacedPatterns(k());
    nthash::SeedNtHash nth(*M_sequence, pattern, 1, k(), M_sequence->length() - k());
    nth.roll();
    for (int i = 0; i < k() - 1; ++i) {
        firstSeedsHashes.push_back(nth.hashes()[i]);
    }
}

const std::string& SequenceInfo::sequence() const {
    return *M_sequence;
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