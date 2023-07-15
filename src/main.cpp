#include <iostream>
#include "SequenceInfo.h"
#include "algorithms.h"
#include "util.h"
#include <chrono>

int main(int argc, char **argv) {

    if (argc != 5) {
        std::cerr << "Wrong argument number." << std::endl;
        return 1;
    }

    size_t k = std::stoi(argv[1]);
    size_t firstK = std::stoi(argv[2]);
    auto reference = util::readFromFASTA(argv[3], true);
    auto r = util::readFromFASTA(argv[4], true);

    auto length = 600000;
    // *r = reference->substr(100000000, length);
    *r = reference->substr(reference->length() - length, length);

    if ((*r)[length - 3] != 'T'){
        (*r)[length - 3] = 'T';
    } else {
        (*r)[length - 3] = 'A';
    }


    std::cout << "Reference length: " << reference->length() << std::endl;
    std::cout << "r length: " << r->length() << std::endl;



    //std::shared_ptr<std::string> reference = std::make_shared<std::string>("RRAACTGTACGGGGGGGGCAAGTGCAAAAAAAAATAAAAAAAA");
    //std::shared_ptr<std::string> r = std::make_shared<std::string>("GGGGGGGGCAAGTGTAAAAAAAAATA");

    //*reference = "GGGAACGTAATTGCGGGTC";
    //*r = "GAACGTAATTGCC";

    //*reference = "ATAACGTAATTGCCGAACGTAATTGCCATACGTAATTGCCAAACGTAATTGCT";
    //*r = "TAATTGCCATACT";

    //std::shared_ptr<std::string> reference = std::make_shared<std::string>("AGCTCGATGACAACTTTGGAC");
    //std::shared_ptr<std::string> r = std::make_shared<std::string>("AAAAAAAAAAAGCTG");


    auto begin = std::chrono::steady_clock::now();
    size_t snpPosition = algorithms::findSNPPosition(reference, r, k, firstK);
    auto end = std::chrono::steady_clock::now();

    if (snpPosition == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at r position " << snpPosition << std::endl;
    }
    std::cout << "Run matching in " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000 << " seconds" << std::endl;
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
