//
// Created by Gabriele Canesi on 01/07/23.
//

#ifndef SPACEDSEEDS_SEQUENCEINFO_H
#define SPACEDSEEDS_SEQUENCEINFO_H

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

class SequenceInfo {
    std::shared_ptr<std::string> M_sequence;
    std::vector<uint64_t> firstSeedsHashes;
    std::unordered_map<uint64_t, std::vector<size_t>> hashTable;
    std::vector<size_t> nullVector;
    size_t M_k;

    SequenceInfo(const std::shared_ptr<std::string> &sequence, const std::vector<std::string> &seeds, bool builhHashTable = true);
    void extendLastSeeds();
public:
    static SequenceInfo buildForReference(const std::shared_ptr<std::string> &sequence, size_t k);
    static SequenceInfo buildForSubstring(const std::shared_ptr<std::string> &sequence, size_t k);

    uint64_t hashAtPosition(size_t position) const;
    const std::vector<size_t>& positionsForHash(uint64_t hash) const;
    const std::string &sequence() const;
    size_t k() const;
};


#endif //SPACEDSEEDS_SEQUENCEINFO_H
