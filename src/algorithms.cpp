//
// Created by Gabriele Canesi on 01/07/23.
//

#include "algorithms.h"
#include <cstddef>
#include "SequenceInfo.h"
#include <iostream>

size_t algorithms::findSNPPosition(const SequenceInfo &referenceInfo, const SequenceInfo &rInfo) {
    size_t k = rInfo.k();

    for (size_t i = 0; i < rInfo.sequence().size(); ++i) {
        for (size_t refPosition : referenceInfo.positionsForHash(rInfo.hashAtPosition(i))) {

            size_t seedStartOffset = 0;
            if (i >= rInfo.sequence().size() - k) {
                seedStartOffset = k - (rInfo.sequence().size() - i);
            }

            if (referenceInfo.sequence()[refPosition + seedStartOffset] != rInfo.sequence()[i]) {

                bool match = true;


                for (size_t j = 1; j * k <= i - seedStartOffset && match; ++j) {
                    if (referenceInfo.hashAtPosition(refPosition - k * j) != rInfo.hashAtPosition(i - k * j - seedStartOffset) ||
                    referenceInfo.sequence()[refPosition - k * j] != rInfo.sequence()[i - k * j - seedStartOffset]) {
                        match = false;
                    }
                }

                for (size_t j = 1; j * k + i < rInfo.sequence().length() - k && match; ++j) {
                    if (referenceInfo.hashAtPosition(refPosition + k * j) != rInfo.hashAtPosition(i + k * j) ||
                    referenceInfo.sequence()[refPosition + k * j] != rInfo.sequence()[i + k * j]) {
                        match = false;
                    }
                }

                if (match && seedStartOffset == 0 && (rInfo.hashAtPosition(rInfo.sequence().length() - k) !=  referenceInfo.hashAtPosition(refPosition - i + (rInfo.sequence().length() - k)) ||
                        rInfo.sequence()[rInfo.sequence().length() - k] != referenceInfo.sequence()[refPosition - i + (rInfo.sequence().length() - k)])) {
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
                    std::cout << "Match: " << refPosition + seedStartOffset << ", " << i << std::endl;
                    //return refPosition + seedStartOffset;
                }
            }
        }
    }

    return -1;
}