//
// Created by Gabriele Canesi on 19/07/23.
//

#include "Test.h"
#include "CliCommand.h"
#include "../algorithms.h"
#include <iostream>

Test::Test(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r, size_t k,
           size_t firstK, double bloomFilterThreshold, double firstThreshold) :
        CliCommand(reference, r, k, firstK, bloomFilterThreshold, firstThreshold) {}



void Test::run() {
    auto length = 6000000;
    *r = reference->substr(reference->length() / 2, length);

    if ((*r)[length - 3] != 'T'){
        (*r)[length - 3] = 'T';
    } else {
        (*r)[length - 3] = 'A';
    }

    auto begin = std::chrono::steady_clock::now();
    size_t position = algorithms::findSNPPosition(reference, r, k, firstK, bloomFilterThreshold, firstThreshold);
    auto end = std::chrono::steady_clock::now();

    if (position == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at r position " << position << std::endl;
    }
    std::cout << "Run matching in " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000 << " seconds" << std::endl;
    assert(position == length - 3);

    *reference = "RRAACTGTACGGGGGGGGCAAGTGCAAAAAAAAATAAAAAAAA";
    *r = "GGGGGGGGCAAGTGTAAAAAAAAATA";
    position = algorithms::findSNPPosition(reference, r, 5, 5, 0.1, 0.7);
    assert(position == 14);

    *reference = "GGGAACGTAATTGCGGGTC";
    *r = "GAACGTAATTGCC";
    *reference = "ATAACGTAATTGCCGAACGTAATTGCCATACGTAATTGCCAAACGTAATTGCT";
    *r = "TAATTGCCATACT";
    *reference = "AGCTCGATGACAACTTTGGAC";
    *r = "AAAAAAAAAAAGCTG";
}