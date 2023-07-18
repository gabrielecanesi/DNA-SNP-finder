//
// Created by Gabriele Canesi on 01/07/23.
//

#include "algorithms.h"
#include "SequenceInfo.h"
#include "SequenceFilter.h"
#include <iostream>

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

    std::cout << similarities[indices[0]] << ", " << indices[0] * size << ", " << indices[0] << ", " << indices.size() << std::endl;
    //std::cout << similarities << std::endl;
    std::cout << "Completed similarities computation" << std::endl;

    auto rInfo = SequenceInfo::buildForSubstring(r->c_str(), k, r->length());
    auto& rKmers = rInfo.exactKmers();
    auto& bloomFilter = rInfo.positionsHashTable();

     for (size_t chunk = 0; chunk < similarities.size() && similarities[indices[chunk]] >= blockThreshold; ++chunk) {

         size_t start = indices[chunk] > 0 ? 1 : 0;
         size_t end = indices[chunk] < indices.size() - 1 ? 1 : 0;
         size_t blockSize = size * (1 + start + end);
         if (indices[chunk] == similarities.size() - 1 || (indices.size() > 1 && indices[chunk] == similarities.size() - 2)) {
             blockSize += reference->length() % size;
         }

         auto referenceInfo = SequenceInfo::buildForReference(reference->c_str() + (indices[chunk] - start) * size, k, blockSize, bloomFilter);
         auto& referenceKmers = referenceInfo.exactKmers();
         std::cout << "Built" << std::endl;


         for (auto& pair : referenceInfo.orderedSubstringHashes()) {
             for (size_t i : rInfo.positionsForHash(pair.first)) {
                 auto &positionsForHash = referenceInfo.positionsForHash(rInfo.hashAtPosition(i));

                 for (size_t refPosition : positionsForHash) {

                     size_t seedStartOffset = 0;
                     if (i >= r->length() - k) {
                         seedStartOffset = k - (r->length() - i);
                     }

                     if ((referenceInfo.sequence())[refPosition + seedStartOffset] !=(*r)[i]) {
                         bool match = true;

                         if (match && seedStartOffset == 0 && refPosition - i + r->length() - k <= referenceKmers.size() && rKmers[r->length() - k] != referenceKmers[refPosition - i + r->length() - k]) {
                             match = false;
                         }

                         if (match && i >= k && refPosition >= (i + seedStartOffset) && rKmers[0] != referenceKmers[refPosition - i + seedStartOffset]) {
                             match = false;
                         }

                         if (i < k) {
                             for (size_t j = 0; match && j < i; ++j) {
                                 if (rInfo.sequence()[j] != referenceInfo.sequence()[refPosition - i + j]) {
                                     match = false;
                                 }
                             }
                         }

                         for (size_t j = 1; match && j * k <= i - seedStartOffset  && k * j <= refPosition; ++j) {
                             if (referenceKmers[refPosition - k * j] != rKmers[i - k * j - seedStartOffset]) {
                                 match = false;
                             }
                         }

                         for (size_t j = 1; match && j * k + i < r->length() - k; ++j) {
                             if (referenceKmers[refPosition + k * j] != rKmers[i + k * j]) {
                                 match = false;
                             }
                         }

                         if (match) {
                             std::cout << "Match: " << (refPosition + seedStartOffset) + size * (indices[chunk] - start) << ", " << i << std::endl;
                             return i;
                         }
                     }
                 }
             }
         }
     }


    return -1;
}





size_t algorithms::findSNPPositionBasic(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r, size_t k) {

    auto rInfo = SequenceInfo::buildForSubstring(r->c_str(), k, r->length());
    auto& rKmers = rInfo.exactKmers();


    auto referenceInfo = SequenceInfo::buildForReference(reference->c_str(), k, reference->length(), rInfo.positionsHashTable());
    auto& referenceKmers = referenceInfo.exactKmers();
    std::cout << "Built" << std::endl;


    for (auto& pair : referenceInfo.orderedSubstringHashes()) {
        for (size_t i : rInfo.positionsForHash(pair.first)) {
            auto &positionsForHash = referenceInfo.positionsForHash(rInfo.hashAtPosition(i));

            for (size_t refPosition : positionsForHash) {

                size_t seedStartOffset = 0;
                if (i >= r->length() - k) {
                    seedStartOffset = k - (r->length() - i);
                }

                if ((referenceInfo.sequence())[refPosition + seedStartOffset] !=(*r)[i]) {
                    bool match = true;

                    if (match && seedStartOffset == 0 && refPosition - i + r->length() - k <= referenceKmers.size() && rKmers[r->length() - k] != referenceKmers[refPosition - i + r->length() - k]) {
                        match = false;
                    }

                    if (match && i >= k && refPosition >= (i + seedStartOffset) && rKmers[0] != referenceKmers[refPosition - i + seedStartOffset]) {
                        match = false;
                    }

                    if (i < k) {
                        for (size_t j = 0; match && j < i; ++j) {
                            if (rInfo.sequence()[j] != referenceInfo.sequence()[refPosition - i + j]) {
                                match = false;
                            }
                        }
                    }

                    for (size_t j = 1; match && j * k <= i - seedStartOffset  && k * j <= refPosition; ++j) {
                        if (referenceKmers[refPosition - k * j] != rKmers[i - k * j - seedStartOffset]) {
                            match = false;
                        }
                    }

                    for (size_t j = 1; match && j * k + i < r->length() - k; ++j) {
                        if (referenceKmers[refPosition + k * j] != rKmers[i + k * j]) {
                            match = false;
                        }
                    }

                    if (match) {
                        std::cout << "Match: " << refPosition + seedStartOffset << ", " << i << std::endl;
                        return i;
                    }
                }
            }
        }
    }

    return -1;
}
