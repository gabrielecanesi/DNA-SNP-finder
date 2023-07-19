//
// Created by Gabriele Canesi on 19/07/23.
//

#ifndef SPACEDSEEDS_FIND_H
#define SPACEDSEEDS_FIND_H


#include "CliCommand.h"

class Find : public CliCommand{

    void run() override;

public:
    Find(const std::shared_ptr<std::string>& reference, const std::shared_ptr<std::string>& r, size_t k, size_t firstK,
         double bloomFilterThreshold, double firstThreshold);

    ~Find() override = default;
};


#endif //SPACEDSEEDS_FIND_H
