//
// Created by Gabriele Canesi on 01/07/23.
//

#include "SequenceInfo.h"
#include "util.h"
#include "nthash/nthash.hpp"

SequenceInfo::SequenceInfo(const char *sequence, const std::vector<std::string> &seeds, size_t length,  const ankerl::unordered_dense::map<uint64_t, std::vector<size_t>>& toCompare, bool buildHashTable) :
        M_sequence(sequence), firstSeedsHashes(), hashTable(), M_k(seeds[0].length()), sequence_length(length) {


    nthash::SeedNtHash nth(sequence, sequence_length, seeds, 1, M_k);

    while (nth.roll()) {
        for (size_t i = 0; i < k(); ++i) {
            if (toCompare.contains(nth.hashes()[i])) {
                util::addToHashTable(nth.hashes()[i], nth.get_pos(), hashTable);
           }
        }
    }

    for (auto& pair : hashTable) {
        orderedRHashes.insert({pair.first, pair.second.size()});
    }

    extractKmers();
}

SequenceInfo::SequenceInfo(const char *sequence, const std::vector<std::string> &seeds, size_t length, bool buildHashTable) :
        M_sequence(sequence), firstSeedsHashes(), hashTable(), M_k(seeds[0].length()), sequence_length(length) {

    nthash::SeedNtHash nth(sequence, sequence_length, seeds, 1, M_k, 0);
    for (size_t i = 0; i <= length - k(); ++i) {
        size_t limit = i < length - k() ? 1 : k();
        nth.roll();
        for (int j = 0; j < limit; ++j) {
            firstSeedsHashes.push_back(nth.hashes()[j]);
            util::addToHashTable(nth.hashes()[j], nth.get_pos() + j, hashTable);
        }
    }
    extractKmers();
}

SequenceInfo SequenceInfo::buildForReference(const char *sequence, size_t k, size_t length, const ankerl::unordered_dense::map<uint64_t, std::vector<size_t>>& rFilter) {
    auto pattern = util::buildAllSpacedPatterns(k);
    return {sequence, pattern, length, rFilter, true};
}

SequenceInfo SequenceInfo::buildForSubstring(const char *sequence, size_t k, size_t length) {
    auto pattern = util::buildAllSpacedPatterns(k);
    return {sequence, pattern, length, false};
}

void SequenceInfo::extractKmers() {
    nthash::NtHash nth(sequence(), sequence_length, 1, k(), 0);
    while (nth.roll()) {
        kmers.push_back(nth.hashes()[0]);
    }
}

const char * SequenceInfo::sequence() const {
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


const std::vector<uint64_t>& SequenceInfo::exactKmers() const {
    return kmers;
}

const std::set<std::pair<uint64_t, size_t>, LessFrequent>& SequenceInfo::orderedSubstringHashes() const {
    return orderedRHashes;
}

const ankerl::unordered_dense::map<uint64_t, std::vector<size_t>>& SequenceInfo::positionsHashTable() const {
    return hashTable;
}