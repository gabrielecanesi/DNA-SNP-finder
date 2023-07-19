//
// Created by Gabriele Canesi on 01/07/23.
//

#ifndef SPACEDSEEDS_ALGORITHMS_H
#define SPACEDSEEDS_ALGORITHMS_H

#include <cstddef>
#include "SequenceInfo.h"

namespace algorithms {
    size_t findSNPPosition(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r,
                           size_t k, size_t firstK, double bloomFilterThreshold, double similarityThreshold);
    size_t findSNPPositionBasic(const std::shared_ptr<std::string> &reference,
                                const std::shared_ptr<std::string> &r,size_t k);
}

#endif //SPACEDSEEDS_ALGORITHMS_H
