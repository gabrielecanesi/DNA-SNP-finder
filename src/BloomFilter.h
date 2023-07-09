//
// Created by Gabriele Canesi on 09/07/23.
//

#ifndef SPACEDSEEDS_BLOOMFILTER_H
#define SPACEDSEEDS_BLOOMFILTER_H


#include <cstddef>
#include <cstdint>
#include <vector>

class BloomFilter {
    bool* array;
    size_t size;
    size_t M_numberOfInserts;
    static std::vector<uint64_t> hashValues;

    size_t hash(uint64_t element, size_t number) const;

public:
    BloomFilter(size_t size);
    ~BloomFilter();
    void insert(uint64_t element);
    bool contains(uint64_t element);
    size_t elementCount() const;
    void clear();
};


#endif //SPACEDSEEDS_BLOOMFILTER_H
