//
// Created by Gabriele Canesi on 01/07/23.
//

#ifndef SPACEDSEEDS_ALGORITHMS_H
#define SPACEDSEEDS_ALGORITHMS_H

#include <cstddef>
#include "SequenceInfo.h"

namespace algorithms {
    size_t findSNPPosition(const SequenceInfo &referenceInfo, const SequenceInfo & rInfo);
}

#endif //SPACEDSEEDS_ALGORITHMS_H
