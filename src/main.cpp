#include <iostream>
#include "SequenceInfo.h"
#include "algorithms.h"
#include "util.h"
#include <chrono>


void test();

int main(int argc, char **argv) {

    if (argc != 6) {
        std::cerr << "Wrong argument count!" << std::endl;
        return 1;
    }

    size_t k = std::stoi(argv[1]);
    size_t firstK = std::stoi(argv[2]);
    auto reference = util::readFromFASTA(argv[3], true);
    auto r = util::readFromFASTA(argv[4], true);
    double bloomFilterThreshold = std::stod(argv[5]);

   auto length = 100000;
    *r = reference->substr(reference->length() / 2, length);

        if ((*r)[length - 3] != 'T'){
        (*r)[length - 3] = 'T';
    } else {
        (*r)[length - 3] = 'A';
    }

    std::cout << "Reference length: " << reference->length() << std::endl;
    std::cout << "r length: " << r->length() << std::endl;

    auto begin = std::chrono::steady_clock::now();
    size_t snpPosition = algorithms::findSNPPosition(reference, r, k, firstK, bloomFilterThreshold);
    //size_t snpPosition = algorithms::findSNPPositionBasic(reference, r, k);
    auto end = std::chrono::steady_clock::now();

    if (snpPosition == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at r position " << snpPosition << std::endl;
    }
    std::cout << "Run matching in " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000 << " seconds" << std::endl;

    return 0;
}


void test() {
    std::shared_ptr<std::string> reference = std::make_shared<std::string>("RRAACTGTACGGGGGGGGCAAGTGCAAAAAAAAATAAAAAAAA");
    std::shared_ptr<std::string> r = std::make_shared<std::string>("GGGGGGGGCAAGTGTAAAAAAAAATA");
    *reference = "GGGAACGTAATTGCGGGTC";
    *r = "GAACGTAATTGCC";
    *reference = "ATAACGTAATTGCCGAACGTAATTGCCATACGTAATTGCCAAACGTAATTGCT";
    *r = "TAATTGCCATACT";
    *reference = "AGCTCGATGACAACTTTGGAC";
    *r = "AAAAAAAAAAAGCTG";
}