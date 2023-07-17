//
// Created by Gabriele Canesi on 09/07/23.
//

#include "BloomFilter.h"
#include <algorithm>
#include <iostream>

//std::vector<uint64_t> BloomFilter::hashValues = {6803, 5347};
std::vector<uint64_t> BloomFilter::hashValues = {31, 53};
BloomFilter::BloomFilter(size_t size) : M_numberOfInserts(0) , size(size) {
    if (size % sizeof(bool) > 0) {
        this->size = ((size + sizeof(bool)) / sizeof(bool)) * sizeof(bool);
    }


    array = new bool[this->size / sizeof(bool)];
    std::fill(array, array + this->size / sizeof(bool), 0);
}

size_t BloomFilter::hash(uint64_t element, size_t number) const {
    auto h1 = (hashValues[0] * element) % (size);
    auto h2 = (number * (hashValues[0] * element)) % (size);
    return (h1 + h2) % size;
}


void BloomFilter::insert(uint64_t element) {
    if (contains(element)) {
        return;
    }

    for (int i = 0; i < 2; ++i) {
        auto position = hash(element, i);
        bool toInsert = 0b10000000 >> (position % sizeof(bool));
        array[position / sizeof(bool)] |= toInsert;
    }
}

bool BloomFilter::contains(uint64_t element) const {
    bool match = true;
    for (size_t i = 0; i < 2 && match; ++i) {
        auto position = hash(element, i);
        bool offset = position % sizeof(bool);
        match = ((array[position / sizeof(bool)] >> (sizeof(bool) - 1 - offset)) & 0b00000001) & match;
    }

    return match;
}

size_t BloomFilter::elementCount() const {
    return M_numberOfInserts;
}

void BloomFilter::clear() {
    std::fill(array, array + size / sizeof(bool), 0);
    M_numberOfInserts = 0;
}


BloomFilter::~BloomFilter() {
    delete[] array;
}

BloomFilter::BloomFilter() : BloomFilter(100){}

