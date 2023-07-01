#include <iostream>
#include "SequenceInfo.h"
#include "algorithms.h"

int main(int argc, char **argv) {
    size_t k = 5;

    
    std::shared_ptr<std::string> reference = std::make_shared<std::string>("RRAACTGTACGGGGGGGGCAAGTGCAAAAAAAAATAAAAAAAA");
    std::shared_ptr<std::string> r = std::make_shared<std::string>("GGGGGGGGCAAGTGTAAAAAAAAATA");

    *reference = "GGGAACGTAATTGCGGGTC";
    *r = "GAACGTAATTGCC";

    *reference = "TAACGTAATTGCCGAACGTAATTGCCATACGTAATTGCCAAACGTAATTGCT";
    *r = "AAACGTAATTGCC";


    auto referenceInfo = SequenceInfo::buildForReference(reference, k);
    auto rInfo = SequenceInfo::buildForSubstring(r, k);

    size_t snpPosition = algorithms::findSNPPosition(referenceInfo, rInfo);

    if (snpPosition == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at reference position " << snpPosition << std::endl;
    }

    return 0;
}
