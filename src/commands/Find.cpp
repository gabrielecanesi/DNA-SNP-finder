//
// Created by Gabriele Canesi on 19/07/23.
//

#include "Find.h"
#include "CliCommand.h"
#include <iostream>
#include "../algorithms.h"

Find::Find(const std::shared_ptr<std::string>& reference, const std::shared_ptr<std::string>& r, size_t k,
           size_t firstK, double bloomFilterThreshold, double firstThreshold) :
           CliCommand(reference, r, k, firstK, bloomFilterThreshold, firstThreshold) {}


void Find::run() {
    auto begin = std::chrono::steady_clock::now();
    size_t snpPosition = algorithms::findSNPPosition(reference, r, k, firstK,
                                                     bloomFilterThreshold, firstThreshold);
    auto end = std::chrono::steady_clock::now();

    if (snpPosition == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at r position " << snpPosition << std::endl;
    }
    std::cout << "Run matching in " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(end -
                            begin).count() / 1000 << " seconds" << std::endl;
}