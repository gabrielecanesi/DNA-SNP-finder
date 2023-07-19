//
// Created by Gabriele Canesi on 01/07/23.
//

#include "algorithms.h"
#include "SequenceInfo.h"
#include "SequenceFilter.h"
#include <iostream>

static inline bool matches(size_t i, size_t k,const SequenceInfo& referenceInfo,const SequenceInfo& rInfo,
                                        const std::vector<uint64_t>& rKmers,
                                        const std::vector<uint64_t>& referenceKmers, size_t size,
                                        size_t start, size_t indexChunk) {
    auto &positionsForHash = referenceInfo
            .positionsForHash(rInfo.hashAtPosition(i));

    for (size_t refPosition : positionsForHash) {


        size_t seedStartOffset = 0;
        if (i >= rInfo.sequenceLength() - k) {
            seedStartOffset = k - (rInfo.sequenceLength() - i);
        }

        if ((referenceInfo.sequence())[refPosition + seedStartOffset] != rInfo.sequence()[i]) {
            bool match = true;

            if (match && seedStartOffset == 0 &&
                refPosition - i + rInfo.sequenceLength() - k <= referenceKmers.size() &&
                rKmers[rInfo.sequenceLength() - k] != referenceKmers[refPosition - i + rInfo.sequenceLength() - k]) {
                match = false;
            }

            if (match && i >= k && refPosition >= (i + seedStartOffset) &&
                rKmers[0] != referenceKmers[refPosition - i + seedStartOffset]) {
                match = false;
            }

            if (i < k) {
                for (size_t j = 0; match && j < i; ++j) {
                    if (rInfo.sequence()[j] != referenceInfo.sequence()[refPosition - i + j]) {
                        match = false;
                    }
                }
            }

            for (size_t j = 1; match && j * k <= i - seedStartOffset && k * j <= refPosition; ++j) {
                if (referenceKmers[refPosition - k * j] != rKmers[i - k * j - seedStartOffset]) {
                    match = false;
                }
            }

            for (size_t j = 1; match && j * k + i < rInfo.sequenceLength() - k; ++j) {
                if (referenceKmers[refPosition + k * j] != rKmers[i + k * j]) {
                    match = false;
                }
            }

            if (match) {
                std::cout << "Match: " << (refPosition + seedStartOffset) +
                                          size * (indexChunk - start) << ", " << i << std::endl;
                return true;
            }
        }
    }
    return false;
}



size_t algorithms::findSNPPosition(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r,
                                   size_t k, size_t firstK, double bloomFilterThreshold, double blockThreshold) {


    size_t size = ceil((double) r->length() / 2);
    SequenceFilter filter(*reference, *r, firstK, size, bloomFilterThreshold);
    auto& similarities = filter.getSimilarities();
    std::vector<size_t> indices(similarities.size(), 0);
    for (int i = 0; i < similarities.size(); ++i) {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return similarities[a] > similarities[b];
    });

    std::cout << similarities[indices[0]] << ", " << indices[0] * size << ", " << indices[0] <<
    ", " << indices.size() << std::endl;

    auto rInfo = SequenceInfo::buildForSubstring(r->c_str(), k, r->length());
    auto& rKmers = rInfo.exactKmers();
    auto& bloomFilter = rInfo.positionsHashTable();

     for (size_t chunk = 0; chunk < similarities.size() && similarities[indices[chunk]] >= blockThreshold; ++chunk) {

         size_t start = indices[chunk] > 0 ? 1 : 0;
         size_t end = indices[chunk] < indices.size() - 1 ? 1 : 0;
         size_t blockSize = size * (1 + start + end);

         if (indices[chunk] == similarities.size() - 1 ||
         (indices.size() > 1 && indices[chunk] == similarities.size() - 2)) {
             blockSize += reference->length() % size;
         }

         auto referenceInfo = SequenceInfo::buildForReference(reference->c_str() +
                 (indices[chunk] - start) * size, k, blockSize, bloomFilter, rKmers);
         auto& referenceKmers = referenceInfo.exactKmers();
         std::cout << "Built block index" << std::endl;


         for (auto& pair : referenceInfo.orderedSubstringHashes()) {
             for (size_t i : rInfo.positionsForHash(pair.first)) {
                 if (matches(i, k, referenceInfo, rInfo, rKmers, referenceKmers,
                                                      size, start, indices[chunk])) {
                     return i;
                 }
             }
         }
     }


    return -1;
}


size_t algorithms::findSNPPositionBasic(const std::shared_ptr<std::string> &reference,
                                        const std::shared_ptr<std::string> &r, size_t k) {

    auto rInfo = SequenceInfo::buildForSubstring(r->c_str(), k, r->length());
    auto& rKmers = rInfo.exactKmers();

    auto referenceInfo = SequenceInfo::buildForReference(reference->c_str(), k,
                                                         reference->length(), rInfo.positionsHashTable(),
                                                         rKmers);
    auto& referenceKmers = referenceInfo.exactKmers();
    std::cout << "Built block index" << std::endl;
    for (auto& pair : referenceInfo.orderedSubstringHashes()) {
        for (size_t i : rInfo.positionsForHash(pair.first)) {
            if (matches(i, k, referenceInfo, rInfo, rKmers, referenceKmers,
                        1, 0, 0)) {
                return i;
            }
        }
    }

    return -1;
}
