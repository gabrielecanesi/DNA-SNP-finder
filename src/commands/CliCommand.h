//
// Created by Gabriele Canesi on 19/07/23.
//

#ifndef SPACEDSEEDS_CLICOMMAND_H
#define SPACEDSEEDS_CLICOMMAND_H
#include <getopt.h>
#include <memory>

class CliCommand {

protected:
    std::shared_ptr<std::string> reference;
    std::shared_ptr<std::string> r;
    size_t k;
    size_t firstK;
    double bloomFilterThreshold;
    double firstThreshold;

    CliCommand(const std::shared_ptr<std::string>& reference,
               const std::shared_ptr<std::string>& r, size_t k, size_t firstK,
               double bloomFilterThreshold, double firstThreshold);

public:
    virtual void run() = 0;
    static CliCommand* create(int argc, char** argv);
    virtual ~CliCommand() = default;
};


#endif //SPACEDSEEDS_CLICOMMAND_H
