//
// Created by Gabriele Canesi on 01/07/23.
//

#include "algorithms.h"
#include <cstddef>
#include "SequenceInfo.h"
#include "SequenceFilter.h"
#include <iostream>
#include "util.h"
#include <nthash/nthash.hpp>

size_t algorithms::findSNPPosition(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r,
                                   size_t k, size_t firstK) {

    size_t size = ceil((double) r->length() / 2);
     SequenceFilter info(*reference, *r, firstK, size);
     auto& similarities = info.getJaccard();
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

     for (size_t chunk = 0; chunk < similarities.size() && similarities[indices[chunk]] >= 0.5; ++chunk) {

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
             //std::cout << pair << std::endl;
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


        auto referenceInfo = SequenceInfo::buildForReference(reference->c_str(), k, reference->length(),{});
        auto rInfo = SequenceInfo::buildForSubstring(r->c_str(), k, r->length());


        for (size_t i = 0; i < r->length(); ++i) {
            auto &positionsForHash = referenceInfo.positionsForHash(rInfo.hashAtPosition(i));
            //std::cout << "Matching positions: " << positionsForHash.size() << std::endl;
            for (size_t refPosition : positionsForHash) {

                size_t seedStartOffset = 0;
                if (i >= r->length() - k) {
                    seedStartOffset = k - (r->length() - i);
                }

                if (referenceInfo.sequence()[refPosition + seedStartOffset] != rInfo.sequence()[i]) {

                    bool match = true;



                    for (size_t j = 1; j * k <= i - seedStartOffset  && k * j <= refPosition && match; ++j) {
                        if (referenceInfo.hashAtPosition(refPosition - k * j) != rInfo.hashAtPosition(i - k * j - seedStartOffset) ||
                            referenceInfo.sequence()[refPosition - k * j] != rInfo.sequence()[i - k * j - seedStartOffset]) {
                            match = false;
                        }
                    }

                    for (size_t j = 1; j * k + i < r->length() - k && match; ++j) {
                        if (referenceInfo.hashAtPosition(refPosition + k * j) != rInfo.hashAtPosition(i + k * j) ||
                            referenceInfo.sequence()[refPosition + k * j] != rInfo.sequence()[i + k * j]) {
                            match = false;
                        }
                    }

                    if (match && seedStartOffset == 0 && (rInfo.hashAtPosition(r->length() - k) !=  referenceInfo.hashAtPosition(refPosition - i + (r->length()- k)) ||
                                                          rInfo.sequence()[r->length() - k] != referenceInfo.sequence()[refPosition - i + (r->length() - k)])) {
                        match = false;
                    }

                    if (match && i >= k && (rInfo.hashAtPosition(0) != referenceInfo.hashAtPosition(refPosition - i + seedStartOffset) ||
                                            rInfo.sequence()[0] != referenceInfo.sequence()[refPosition - i + seedStartOffset])) {
                        match = false;
                    }

                    if (i < k) {
                        for (size_t j = 0; j < i && match; ++j) {
                            if (rInfo.sequence()[j] != referenceInfo.sequence()[refPosition - i + j]) {
                                match = false;
                            }
                        }
                    }

                    if (match) {
                        std::cout << "Match: " << (refPosition + seedStartOffset - r->length()) << ", " << i << std::endl;
                        //return refPosition + seedStartOffset;
                    }
                }
            }
        }

    return -1;
}



size_t algorithms::findSNPPositionApprox(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r, size_t k) {

    auto size = r->length();
    SequenceFilter info(*reference, *r, k, size);
    auto jaccardsWithR = info.getJaccard();
    std::vector<size_t> indices(jaccardsWithR.size(), 0);
    for (int i = 1; i < jaccardsWithR.size(); ++i) {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return jaccardsWithR[a] > jaccardsWithR[b];
    });


    std::cout << jaccardsWithR[indices[0]] << ", " << indices[0] * r->length() << std::endl;
    std::cout << jaccardsWithR << std::endl;

    std::cout << "Completed similarities computation" << std::endl;

    auto spacedSeeds = util::buildAllSpacedPatterns(k);
    auto seedHashes = new uint64_t[(r->length() - k + 1) * k];
    auto rKmers = info.rKmers();
    nthash::SeedNtHash rNtHash(*r, spacedSeeds, 1, k);

    for (size_t j = 0; j < r->length() - k + 1; ++j) {
        rNtHash.roll();
        for (size_t i = 0; i < k; ++i) {
            seedHashes[i * (r->length() - k + 1) + j] = rNtHash.hashes()[i];
        }
    }
    std::cout << "Built r spaced kmers" << std::endl;

    BloomFilter referenceHashes(1000000000);

    for (size_t chunk = 0; chunk < jaccardsWithR.size() && jaccardsWithR[indices[chunk]] >= 0.4; ++chunk) {
        referenceHashes.clear();
        size_t start = indices[chunk] > 0 ? 1 : 0;
        size_t end = indices[chunk] < indices.size() - 1 ? 1 : 0;

        nthash::NtHash ntHashKmerReference(reference->c_str() + (indices[chunk] - start) * size, (1 + start + end) * size, 1, k, 0);

        while (ntHashKmerReference.roll()) {
            referenceHashes.insert(ntHashKmerReference.hashes()[0]);
        }

        nthash::SeedNtHash nthSeedReference(reference->c_str() + (indices[chunk] - start) * size, (1 + start + end) * size, spacedSeeds, 1, k, 0);
        while (nthSeedReference.roll()) {
            for (size_t i = 0; i < k; ++i) {
                referenceHashes.insert(nthSeedReference.hashes()[i]);
            }
        }
        std::cout << "Built" << std::endl;

        for (size_t mutIndex = 0; mutIndex < r->length() - k + 1; ++mutIndex) {
            if (!referenceHashes.contains(rKmers[mutIndex])) {
                return mutIndex + k - 1;
            }
        }

    }

    return -1;
}
