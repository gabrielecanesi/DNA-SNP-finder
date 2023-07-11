#include <iostream>
#include "SequenceInfo.h"
#include "algorithms.h"
#include "util.h"
#include <chrono>

int main(int argc, char **argv) {

    if (argc != 4) {
        throw std::runtime_error("Wrong argument number");
    }

    size_t k = std::stoi(argv[1]);
    auto reference = util::readFromFASTA(argv[2], true);
    auto r = util::readFromFASTA(argv[3], true);

    auto length = 6000000;
    *r = reference->substr(100000000, length);

    if ((*r)[length - 3] != 'T'){
        (*r)[length - 3] = 'T';
    } else {
        (*r)[length - 3] = 'A';
    }


    std::cout << "Reference length: " << reference->length() << std::endl;

    //std::shared_ptr<std::string> reference = std::make_shared<std::string>("RRAACTGTACGGGGGGGGCAAGTGCAAAAAAAAATAAAAAAAA");
    //std::shared_ptr<std::string> r = std::make_shared<std::string>("GGGGGGGGCAAGTGTAAAAAAAAATA");

   /* *reference = "GGGAACGTAATTGCGGGTC";
    *r = "GAACGTAATTGCC";

    *reference = "ATAACGTAATTGCCGAACGTAATTGCCATACGTAATTGCCAAACGTAATTGCT";
    *r = "TAATTGCCATACT";*/


    auto begin = std::chrono::steady_clock::now();
    size_t snpPosition = algorithms::findSNPPosition(reference, r, k);
    auto end = std::chrono::steady_clock::now();

    if (snpPosition == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at r position " << snpPosition << std::endl;
    }
    std::cout << "Run matching in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " seconds" << std::endl;
    /*begin = std::chrono::steady_clock::now();
    snpPosition = algorithms::findSNPPositionBasic(reference, r, k);
    if (snpPosition == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at reference position " << snpPosition << std::endl;
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Run matching in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " seconds" << std::endl;
*/
    return 0;
}
