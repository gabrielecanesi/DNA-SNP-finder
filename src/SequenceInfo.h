//
// Created by Gabriele Canesi on 01/07/23.
//

#ifndef SPACEDSEEDS_SEQUENCEINFO_H
#define SPACEDSEEDS_SEQUENCEINFO_H

#include <memory>
#include <string>
#include <vector>
#include <list>
#include <ankerl/unordered_dense.h>

class SequenceInfo {
    const char *M_sequence;
    size_t sequence_length;
    std::vector<uint64_t> firstSeedsHashes;
    ankerl::unordered_dense::map<uint64_t, std::list<size_t>> hashTable;
    std::list<size_t> nullVector;
    size_t M_k;

    SequenceInfo(const char *sequence, const std::vector<std::string> &seeds, size_t length, bool buildHashTable = true);
    void extendLastSeeds();
public:
    static SequenceInfo buildForReference(const char *sequence, size_t k, size_t length);
    static SequenceInfo buildForSubstring(const char *sequence, size_t k, size_t length);

    uint64_t hashAtPosition(size_t position) const;
    const std::list<size_t>& positionsForHash(uint64_t hash) const;
    const char *sequence();
    size_t k() const;
    size_t sequenceLength() const;
};


#endif //SPACEDSEEDS_SEQUENCEINFO_H
