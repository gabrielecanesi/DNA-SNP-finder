#include <iostream>
#include "SequenceInfo.h"
#include "algorithms.h"
#include "util.h"
#include <chrono>


void test(const std::shared_ptr<std::string>& reference, const std::shared_ptr<std::string>& r, size_t k,
          size_t firstK, double bloomFilterThreshold, double blockThreshold);

int main(int argc, char **argv) {

    if (argc != 8) {
        std::cerr << "Wrong argument count!" << std::endl;
        return 1;
    }

    size_t k = std::stoi(argv[2]);
    size_t firstK = std::stoi(argv[3]);
    auto reference = util::readFromFASTA(argv[4], true);
    auto r = util::readFromFASTA(argv[5], true);
    double bloomFilterThreshold = std::stod(argv[6]);
    double similarityThreshold = std::stod(argv[7]);

    std::string cmd = argv[1];

    if (cmd == "test") {
        test(reference, r, k, firstK, bloomFilterThreshold, similarityThreshold);
        return 0;
    }

    auto begin = std::chrono::steady_clock::now();
    size_t snpPosition = algorithms::findSNPPosition(reference, r, k, firstK, bloomFilterThreshold, similarityThreshold);
    auto end = std::chrono::steady_clock::now();


    if (snpPosition == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at r position " << snpPosition << std::endl;
    }
    std::cout << "Run matching in " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000 << " seconds" << std::endl;


    return 0;
}


void test(const std::shared_ptr<std::string>& reference, const std::shared_ptr<std::string>& r, size_t k,
          size_t firstK, double bloomFilterThreshold, double blockThreshold) {


    auto length = 60000000;
    *r = reference->substr(reference->length() / 2, length);

    if ((*r)[length - 3] != 'T'){
        (*r)[length - 3] = 'T';
    } else {
        (*r)[length - 3] = 'A';
    }

    auto begin = std::chrono::steady_clock::now();
    size_t position = algorithms::findSNPPosition(reference, r, k, firstK, bloomFilterThreshold, blockThreshold);
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