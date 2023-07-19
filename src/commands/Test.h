//
// Created by Gabriele Canesi on 19/07/23.
//

#ifndef SPACEDSEEDS_TEST_H
#define SPACEDSEEDS_TEST_H


#include "CliCommand.h"

class Test : public CliCommand {

    void run() override;

public:
    Test(const std::shared_ptr<std::string>& reference, const std::shared_ptr<std::string>& r, size_t k, size_t firstK,
         double bloomFilterThreshold, double firstThreshold);

    ~Test() override = default;
};


#endif //SPACEDSEEDS_TEST_H
