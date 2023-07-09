//
// Created by Gabriele Canesi on 01/07/23.
//

#include "algorithms.h"
#include <cstddef>
#include "SequenceInfo.h"
#include "MinHashInfo.h"
#include <iostream>

template <typename  T>
std::ostream& operator<<(std::ostream &out, const std::vector<T> &vector) {
    out << "[";
    for (size_t i = 0; i < vector.size() - 1; ++i) {
        out << vector[i] << ", ";
    }
    if (vector.size() > 0) {
        out << vector[vector.size() - 1];
    }
    out << "]";

    return out;
}


size_t algorithms::findSNPPosition(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r, size_t k) {

     MinHashInfo info(*reference, *r, k);
     info.computeSignature();
     auto jaccardsWithR = info.computeJaccard();
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

    auto rInfo = SequenceInfo::buildForSubstring(r->c_str(), k, r->length());
     for (size_t chunk = 0; chunk < jaccardsWithR.size() && jaccardsWithR[indices[chunk]] >= 0.4; ++chunk) {

         size_t start = indices[chunk] > 0 ? 1 : 0;
         size_t end = indices[chunk] < indices.size() - 1 ? 1 : 0;

         auto referenceInfo = SequenceInfo::buildForReference(reference->c_str() + (indices[chunk] - start) * r->length(), k, r->length() * (1 + start + end));
         std::cout << "Built" << std::endl;


         for (size_t i = 0; i < rInfo.sequenceLength(); ++i) {
             auto &positionsForHash = referenceInfo.positionsForHash(rInfo.hashAtPosition(i));
             //std::cout << "Matching positions: " << positionsForHash.size() << std::endl;
             for (size_t refPosition : positionsForHash) {

                 size_t seedStartOffset = 0;
                 if (i >= rInfo.sequenceLength() - k) {
                     seedStartOffset = k - (rInfo.sequenceLength() - i);
                 }

                 if (referenceInfo.sequence()[refPosition + seedStartOffset] != rInfo.sequence()[i]) {

                     bool match = true;


                     for (size_t j = 1; j * k <= i - seedStartOffset  && k * j <= refPosition && match; ++j) {
                         if (referenceInfo.hashAtPosition(refPosition - k * j) != rInfo.hashAtPosition(i - k * j - seedStartOffset) ||
                             referenceInfo.sequence()[refPosition - k * j] != rInfo.sequence()[i - k * j - seedStartOffset]) {
                             match = false;
                         }
                     }

                     for (size_t j = 1; j * k + i < rInfo.sequenceLength() - k && match; ++j) {
                         if (referenceInfo.hashAtPosition(refPosition + k * j) != rInfo.hashAtPosition(i + k * j) ||
                             referenceInfo.sequence()[refPosition + k * j] != rInfo.sequence()[i + k * j]) {
                             match = false;
                         }
                     }

                     if (match && seedStartOffset == 0 && (rInfo.hashAtPosition(rInfo.sequenceLength() - k) !=  referenceInfo.hashAtPosition(refPosition - i + (rInfo.sequenceLength()- k)) ||
                                                           rInfo.sequence()[rInfo.sequenceLength() - k] != referenceInfo.sequence()[refPosition - i + (rInfo.sequenceLength() - k)])) {
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
                         std::cout << "Match: " << (refPosition + seedStartOffset) + r->length() * indices[chunk] << ", " << i << std::endl;
                         return i;
                     }
                 }
             }
         }

     }


    return -1;
}





size_t algorithms::findSNPPositionBasic(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r, size_t k) {


        auto referenceInfo = SequenceInfo::buildForReference(reference->c_str(), k, reference->length());
        auto rInfo = SequenceInfo::buildForSubstring(r->c_str(), k, r->length());


        for (size_t i = 0; i < rInfo.sequenceLength(); ++i) {
            auto &positionsForHash = referenceInfo.positionsForHash(rInfo.hashAtPosition(i));
            //std::cout << "Matching positions: " << positionsForHash.size() << std::endl;
            for (size_t refPosition : positionsForHash) {

                size_t seedStartOffset = 0;
                if (i >= rInfo.sequenceLength() - k) {
                    seedStartOffset = k - (rInfo.sequenceLength() - i);
                }

                if (referenceInfo.sequence()[refPosition + seedStartOffset] != rInfo.sequence()[i]) {

                    bool match = true;


                    for (size_t j = 1; j * k <= i - seedStartOffset  && k * j <= refPosition && match; ++j) {
                        if (referenceInfo.hashAtPosition(refPosition - k * j) != rInfo.hashAtPosition(i - k * j - seedStartOffset) ||
                            referenceInfo.sequence()[refPosition - k * j] != rInfo.sequence()[i - k * j - seedStartOffset]) {
                            match = false;
                        }
                    }

                    for (size_t j = 1; j * k + i < rInfo.sequenceLength() - k && match; ++j) {
                        if (referenceInfo.hashAtPosition(refPosition + k * j) != rInfo.hashAtPosition(i + k * j) ||
                            referenceInfo.sequence()[refPosition + k * j] != rInfo.sequence()[i + k * j]) {
                            match = false;
                        }
                    }

                    if (match && seedStartOffset == 0 && (rInfo.hashAtPosition(rInfo.sequenceLength() - k) !=  referenceInfo.hashAtPosition(refPosition - i + (rInfo.sequenceLength()- k)) ||
                                                          rInfo.sequence()[rInfo.sequenceLength() - k] != referenceInfo.sequence()[refPosition - i + (rInfo.sequenceLength() - k)])) {
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
                        std::cout << "Match: " << (refPosition + seedStartOffset) << ", " << i << std::endl;
                        //return refPosition + seedStartOffset;
                    }
                }
            }
        }

    return -1;
}



size_t algorithms::findSNPPositionApprox(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r, size_t k) {

    MinHashInfo info(*reference, *r, k);
    info.computeSignature();
    auto jaccardsWithR = info.computeJaccard();
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

    auto rInfo = SequenceInfo::buildForSubstring(r->c_str(), k, r->length());
    for (size_t chunk = 0; chunk < jaccardsWithR.size() && jaccardsWithR[indices[chunk]] >= 0.4; ++chunk) {

        size_t start = indices[chunk] > 0 ? 1 : 0;
        size_t end = indices[chunk] < indices.size() - 1 ? 1 : 0;

        auto referenceInfo = SequenceInfo::buildForReference(reference->c_str() + (indices[chunk] - start) * r->length(), k, r->length() * (1 + start + end));
        std::cout << "Built" << std::endl;


        for (size_t i = 0; i < rInfo.sequenceLength(); ++i) {
            auto &positionsForHash = referenceInfo.positionsForHash(rInfo.hashAtPosition(i));
            //std::cout << "Matching positions: " << positionsForHash.size() << std::endl;
            for (size_t refPosition : positionsForHash) {

                size_t seedStartOffset = 0;
                if (i >= rInfo.sequenceLength() - k) {
                    seedStartOffset = k - (rInfo.sequenceLength() - i);
                }

                if (referenceInfo.sequence()[refPosition + seedStartOffset] != rInfo.sequence()[i]) {

                    bool match = true;


                    for (size_t j = 1; j * k <= i - seedStartOffset  && k * j <= refPosition && match; ++j) {
                        if (referenceInfo.hashAtPosition(refPosition - k * j) != rInfo.hashAtPosition(i - k * j - seedStartOffset) ||
                            referenceInfo.sequence()[refPosition - k * j] != rInfo.sequence()[i - k * j - seedStartOffset]) {
                            match = false;
                        }
                    }

                    for (size_t j = 1; j * k + i < rInfo.sequenceLength() - k && match; ++j) {
                        if (referenceInfo.hashAtPosition(refPosition + k * j) != rInfo.hashAtPosition(i + k * j) ||
                            referenceInfo.sequence()[refPosition + k * j] != rInfo.sequence()[i + k * j]) {
                            match = false;
                        }
                    }

                    if (match && seedStartOffset == 0 && (rInfo.hashAtPosition(rInfo.sequenceLength() - k) !=  referenceInfo.hashAtPosition(refPosition - i + (rInfo.sequenceLength()- k)) ||
                                                          rInfo.sequence()[rInfo.sequenceLength() - k] != referenceInfo.sequence()[refPosition - i + (rInfo.sequenceLength() - k)])) {
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
                        std::cout << "Match: " << (refPosition + seedStartOffset) + r->length() * indices[chunk] << ", " << i << std::endl;
                        return i;
                    }
                }
            }
        }

    }


    return -1;
}
