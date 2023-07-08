//
// Created by Gabriele Canesi on 01/07/23.
//


#include "util.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>


std::vector<std::string> util::buildAllSpacedPatterns(size_t k) {
    std::vector<std::string> result(k, std::string(k, '1'));
    for (size_t i = 0; i < k; ++i) {
        result[i][i] = '0';
    }

    return result;
}


std::vector<std::string> util::buildFirstSpacedPattern(size_t k) {
    std::vector<std::string> result(1, std::string(k, '1'));
    result[0][0] = '0';

    return result;
}

std::vector<std::string> util::buildRemainingSpacedPatterns(size_t k) {
    std::vector<std::string> result(k - 1, std::string(k, '1'));
    for (size_t i = 1; i < k; ++i) {
        result[i - 1][i] = '0';
    }

    return result;
}



std::shared_ptr<std::string> util::readFromFASTA(const std::string &path, bool skipNonACGT) {
    std::ifstream fileStream(path);
    std::shared_ptr<std::string> line = std::make_shared<std::string>();
    std::getline(fileStream, *line);
    if (line->length() == 0 || (*line)[0] != '>') {
        throw std::runtime_error("Wrong file format");
    }

    *line = "";
    char c;

    while (fileStream.get(c) && c != '>') {
        if ((!skipNonACGT || c == 'A' || c == 'C' || c == 'G' || c == 'T') && c != '\n' && c != '\0') {
            line->push_back(c);
        }
    }

    return line;
}