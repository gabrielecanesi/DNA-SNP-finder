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
#include <set>

struct LessFrequent {
    bool operator()(const std::pair<uint64_t, size_t>& a, const std::pair<uint64_t, size_t>& b) const {
        return a.second <= b.second;
    }
};

class SequenceInfo {
    const char *M_sequence;
    size_t sequence_length;
    std::vector<uint64_t> firstSeedsHashes;
    ankerl::unordered_dense::map<uint64_t, std::vector<size_t>> hashTable;
    std::vector<size_t> nullVector;
    size_t M_k;
    std::set<std::pair<uint64_t, size_t>, LessFrequent> orderedRHashes;
    std::vector<uint64_t> kmers;


    SequenceInfo(const char *sequence, const std::vector<std::string> &seeds, size_t length,
                 const ankerl::unordered_dense::map<uint64_t, std::vector<size_t>>& toCompare,
                 const std::vector<uint64_t>& rKmers, bool buildHashTable = true);
    SequenceInfo(const char *sequence, const std::vector<std::string> &seeds, size_t length, bool buildHashTable = true);
    void extractKmers();
public:
    static SequenceInfo buildForReference(const char *sequence, size_t k, size_t length,
                                          const ankerl::unordered_dense::map<uint64_t, std::vector<size_t>>& rFilter,
                                          const std::vector<uint64_t>& rKmers);
    static SequenceInfo buildForSubstring(const char *sequence, size_t k, size_t length);

    uint64_t hashAtPosition(size_t position) const;
    const std::vector<size_t>& positionsForHash(uint64_t hash) const;
    const char *sequence() const;
    size_t k() const;
    size_t sequenceLength() const;
    const std::vector<uint64_t>& exactKmers() const;
    const std::set<std::pair<uint64_t, size_t>, LessFrequent>& orderedSubstringHashes() const;
    const ankerl::unordered_dense::map<uint64_t, std::vector<size_t>>& positionsHashTable() const;
};


#endif //SPACEDSEEDS_SEQUENCEINFO_H
