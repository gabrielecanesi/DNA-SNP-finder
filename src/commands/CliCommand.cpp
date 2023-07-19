//
// Created by Gabriele Canesi on 19/07/23.
//

#include "CliCommand.h"
#include "../util.h"
#include "Find.h"
#include "Test.h"

CliCommand *CliCommand::create(int argc, char **argv) {
    option commandOptions[4] = {
            {"k", required_argument, nullptr, 'k'},
            {"firstK", required_argument, nullptr, 'f'},
            {"bloomFilterThreshold", required_argument, nullptr, 'b'},
            {"firstThreshold", required_argument, nullptr, 't'}
    };

    if (argc < 3) {
        throw std::invalid_argument("Invalid argument number.");
    }
    auto reference = util::readFromFASTA(argv[2], true);
    auto r = util::readFromFASTA(argv[3], true);
    size_t k = 10;
    size_t firstK = 15;
    double bloomFilterThreshold = 0.1;
    double firstThreshold = 0.8;

    std::string command = argv[1];
    int option;
    while ((option = getopt_long(argc, argv, "k:f:b:t", commandOptions, nullptr)) != -1) {
        auto value = std::string(optarg);
        value.erase(std::remove(value.begin(), value.end(), '='));
        value.erase(std::remove(value.begin(), value.end(), ' '));

        switch (option) {
            case 'k':
                k = std::stoi(value);
                break;
            case 'f':
                firstK = std::stoi(value);
                break;
            case 'b':
                bloomFilterThreshold = std::stod(value);
                break;
            case 't':
                firstThreshold = std::stod(value);
        }
    }

    if (command == "run") {
        return new Find(reference, r, k, firstK, bloomFilterThreshold, firstThreshold);
    } else if (command == "test") {
        return new Test(reference, r, k, firstK, bloomFilterThreshold, firstThreshold);
    }

    throw std::invalid_argument("Command not recognized!");
}


CliCommand::CliCommand(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r, size_t k,
                       size_t firstK, double bloomFilterThreshold, double firstThreshold) :
                       reference(reference), r(r), k(k), firstK(firstK), bloomFilterThreshold(bloomFilterThreshold),
                       firstThreshold(firstThreshold){}
