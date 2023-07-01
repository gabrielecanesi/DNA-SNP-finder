#include <iostream>
#include <nthash/nthash.hpp>


std::vector<std::string> buildFirstAndKmerPatterns(int k) {
    std::vector<std::string> result(2, std::string(k, '1'));
    result[0][0] = '0';
    return result;
}

std::vector<std::string> buildIntermediateSpacedPatterns(int k) {
    std::vector<std::string> result(k - 1, std::string(k, '1'));
    for (int i = 0; i < k - 1; ++i) {
        result[i][i + 1] = '0';
    }

    return result;
}

std::vector<std::string> buildAllPatterns(int k) {
    std::vector<std::string> seeds;
    for (int i = 0; i < k; ++i) {
        seeds.emplace_back(k, '1');
        seeds[i][i] = '0';
    }
    seeds.emplace_back(k, '1');
    return seeds;
}

void addToHashTable(uint64_t hash, size_t position, std::unordered_map<uint64_t, std::vector<size_t>> &hashTable) {
    auto iterator = hashTable.find(hash);
    if (iterator == hashTable.end()) {
        hashTable[hash] = std::vector<size_t>(1, position);
    } else {
        hashTable[hash].push_back(position);
    }
}

struct SequenceInfo {
    size_t k;
    size_t sequenceLength;
    std::unordered_map<uint64_t, std::vector<size_t>> hashTable;
    std::vector<std::vector<uint64_t>> seedsArrays;
    const std::string &sequence;

    SequenceInfo(const std::string &sequence, int k, bool reference) : sequence(sequence), sequenceLength(sequence.length()), k(k),
                                                seedsArrays(std::vector<std::vector<uint64_t>>(k + 1, {})) {

        if (reference) {
            buildReference();
        } else {
            buildSubstring();
        }
    }

private:
    void buildReference() {
        auto seeds = buildAllPatterns(k);
        nthash::SeedNtHash nth(sequence, seeds, 1, k);
        while(nth.roll()) {
            for (int i = 0; i <= k; ++i) {
                addToHashTable(nth.hashes()[i], nth.get_pos(), hashTable);
                seedsArrays[i].push_back(nth.hashes()[i]);
            }
        }
    }

    void buildSubstring() {
        auto intermediateSeeds = buildIntermediateSpacedPatterns(k);
        auto firstAndExactKmers = buildFirstAndKmerPatterns(k);
        nthash::SeedNtHash nth(sequence, firstAndExactKmers, 1, k);
        auto nthIntermediate = nthash::SeedNtHash(sequence, intermediateSeeds, 1, k, sequence.length() - k);

        while(nth.roll()) {
            for (int i = 0; i < 2; ++i) {
                addToHashTable(nth.hashes()[i], nth.get_pos(), hashTable);
                seedsArrays[i * k].push_back(nth.hashes()[i]);
            }
        }

        while(nthIntermediate.roll()) {
            for (int i = 0; i < k - 1; ++i) {
                addToHashTable(nthIntermediate.hashes()[i], nthIntermediate.get_pos(), hashTable);
                seedsArrays[i + 1].push_back(nthIntermediate.hashes()[i]);
            }
        }
    }
};

size_t findSNP(SequenceInfo &referenceInfo, SequenceInfo &rInfo) {
    for (long i = 0; i < rInfo.sequenceLength - rInfo.k + 1; ++i) {

        auto referencePositions = referenceInfo.hashTable.find(rInfo.seedsArrays[0][i]);
        
        if (referencePositions != referenceInfo.hashTable.end()) {
            for (auto position : referencePositions->second) {
                if (referenceInfo.sequence[position] != rInfo.sequence[i]) { // k-mer differs only by the selected position i on r
                    bool match = true;
                    for (long j = 1; i - j * (long) rInfo.k >= 0 && match; ++j) { // Compare all the preceeding kmer hashes

                        if (rInfo.seedsArrays[rInfo.k][i - j * (long) rInfo.k] !=referenceInfo.seedsArrays[referenceInfo.k][position - j * (long) rInfo.k]) {
                            match = false;
                        }
                    }

                    for (long j = 0; j < position % rInfo.k && j < i && match; ++j) {
                        if (rInfo.sequence[j] != referenceInfo.sequence[position - i + j]) {
                            match = false;
                        }
                    }

                    for (int j = 1; i + j * (long) rInfo.k <= rInfo.sequenceLength - rInfo.k && match; ++j) {
                        if (rInfo.seedsArrays[rInfo.k][i + j * (long) rInfo.k] != referenceInfo.seedsArrays[referenceInfo.k][position + j * (long) rInfo.k]) {
                            match = false;
                        }
                    }

                    size_t leftSymbols = (rInfo.sequenceLength - i) % rInfo.k;

                    for (int j = 1; j <= leftSymbols && match; ++j) {
                        if (rInfo.sequence[rInfo.sequenceLength - j] != referenceInfo.sequence[position + (rInfo.sequenceLength - i) - j]) {
                            match = false;
                        }
                    }

                    if (match) {
                        //return position;
                        std::cout << position << " " << i << std::endl;
                    }
                }
            }
        }
    }


    for (int i = 1; i < rInfo.k; ++i) {
        size_t startPosition = rInfo.sequenceLength - rInfo.k;
        size_t hash = rInfo.seedsArrays[i][0];
        auto positions = referenceInfo.hashTable.find(hash);

        if (positions != referenceInfo.hashTable.end()) {
            for (auto position : positions->second) {
                if (referenceInfo.sequence[position + i] != rInfo.sequence[startPosition + i]) {
                    bool match = true;
                    for (long j = 1; (long) startPosition - j * (long) rInfo.k >= 0 && match; ++j) { // Compare all the preceeding kmer hashes

                        if (rInfo.seedsArrays[rInfo.k][startPosition - j * (long) rInfo.k] != referenceInfo.seedsArrays[referenceInfo.k][position - (long) j * rInfo.k]) {
                            match = false;
                        }
                    }

                    for (long j = 0; j < startPosition % rInfo.k && match; ++j) {
                        if (rInfo.sequence[j] != referenceInfo.sequence[position - startPosition + j]) {
                            match = false;
                        }
                    }

                    if (match) {
                        //return position + i;
                        std::cout << position + i << " - " << startPosition + i << std::endl;
                    }
                }
            }
        }
    }

    return -1;
}




int main(int argc, char **argv) {
     int k = 5;
    
    std::string reference = "RRAACTGTACGGGGGGGGCAAGTGCAAAAAAAAATAAAAAAAA";
    std::string r = "GGGGGGGGCAAGTGTAAAAAAAAATA";

    reference = "GGGAACGTAATTGCGGGTC";
    r = "GAACGTAATTGCC";


    SequenceInfo referenceInfo(reference, k, true);
    SequenceInfo rInfo(r, k, false);

    size_t snpPosition = findSNP(referenceInfo, rInfo);

    if (snpPosition == -1) {
        std::cout << "SNP not found." << std::endl;
    } else {
        std::cout << "Found SNP at reference position " << snpPosition << std::endl;
    }

    return 0;
}
