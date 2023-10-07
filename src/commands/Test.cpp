//
// Created by Gabriele Canesi on 19/07/23.
//

#include "Test.h"
#include "CliCommand.h"
#include "../algorithms.h"
#include <iostream>
#include <fstream>
#include <cassert>

Test::Test(const std::shared_ptr<std::string> &reference, const std::shared_ptr<std::string> &r, size_t k,
           size_t firstK, double bloomFilterThreshold, double firstThreshold) :
        CliCommand(reference, r, k, firstK, bloomFilterThreshold, firstThreshold) {}



void Test::run() {
    auto originalReference = std::make_shared<std::string>(*reference);
    auto originalRead = std::make_shared<std::string>(*r);
    auto readLengths = {1000, 10000, 50000, 100000, 500000, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000};
    auto referenceLengths = {100000, 500000, 1000000, 5000000, 10000000, 20000000, 30000000, 40000000,
                             50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 120000000, 140000000, 160000000, 180000000, 200000000, (int) reference->length()};

    std::ofstream file("timingsR.csv");
    file << "read,reference,time\n";
    for (auto lR : referenceLengths) {
        auto l = 100000;
        if (l < lR) {

            auto pos = rand() % l;
            std::cout << "Pos: " << pos << std::endl;
            auto startPos = rand() % (lR - l);
            *r = reference->substr(startPos, l);

            if ((*r)[pos] != 'T'){
                (*r)[pos] = 'T';
            } else {
                (*r)[pos] = 'A';
            }

            std::cout << "Reference length: " << lR << std::endl;
            std::cout << "r length: " << l << std::endl;


            auto begin = std::chrono::steady_clock::now();
            size_t position = algorithms::findSNPPosition(reference, r, k, firstK, bloomFilterThreshold,
                                                          firstThreshold);
            auto end = std::chrono::steady_clock::now();

            if (position == -1) {
                std::cout << "SNP not found." << std::endl;
            } else {
                std::cout << "Found SNP at r position " << position << std::endl;
            }
            assert(position == pos);
            std::cout << "Run matching in " <<
            (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000 <<
            " seconds" << std::endl;

            file << l << "," << lR << "," << (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000 << "\n";
            *reference = *originalReference;
            *r = *originalRead;
        }
    }

    auto newFile = std::ofstream("timingsRead.csv");
    newFile << "read,reference,time\n";
    for (auto l : readLengths) {
        auto lR = reference->length();
        if (l < lR) {

            auto pos = rand() % l;
            std::cout << "Pos: " << pos << std::endl;
            auto startPos = rand() % (lR - l);
            *r = reference->substr(startPos, l);

            if ((*r)[pos] != 'T'){
                (*r)[pos] = 'T';
            } else {
                (*r)[pos] = 'A';
            }

            std::cout << "Reference length: " << lR << std::endl;
            std::cout << "r length: " << l << std::endl;


            auto begin = std::chrono::steady_clock::now();
            size_t position = algorithms::findSNPPosition(reference, r, k, firstK, bloomFilterThreshold,
                                                          firstThreshold);
            auto end = std::chrono::steady_clock::now();

            if (position == -1) {
                std::cout << "SNP not found." << std::endl;
            } else {
                std::cout << "Found SNP at r position " << position << std::endl;
            }
            assert(position == pos);
            std::cout << "Run matching in " <<
                      (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000 <<
                      " seconds" << std::endl;

            newFile << l << "," << lR << "," << (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000 << "\n";

        }
    }


    assert(position == length - 3);

    *reference = "RRAACTGTACGGGGGGGGCAAGTGCAAAAAAAAATAAAAAAAA";
    *r = "GGGGGGGGCAAGTGTAAAAAAAAATA";
    size_t position = algorithms::findSNPPosition(reference, r, 5, 5, 0.1, 0.7);
    assert(position == 14);

    *reference = "GGGAACGTAATTGCGGGTC";
    *r = "GAACGTAATTGCC";
    *reference = "ATAACGTAATTGCCGAACGTAATTGCCATACGTAATTGCCAAACGTAATTGCT";
    *r = "TAATTGCCATACT";
    *reference = "AGCTCGATGACAACTTTGGAC";
    *r = "AAAAAAAAAAAGCTG";
}