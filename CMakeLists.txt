CMAKE_MINIMUM_REQUIRED(VERSION 3.25)

SET (CMAKE_CXX_STANDARD 17)

PROJECT(spacedSeeds)

find_package(unordered_dense CONFIG REQUIRED)
ADD_EXECUTABLE(spacedSeeds src/main.cpp src/SequenceInfo.cpp src/SequenceInfo.h src/util.h src/util.cpp src/algorithms.h src/algorithms.cpp src/SequenceFilter.cpp src/SequenceFilter.h bloom/bloom_filter.hpp src/commands/CliCommand.cpp src/commands/CliCommand.h src/commands/Test.cpp src/commands/Test.h src/commands/Find.cpp src/commands/Find.h)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/ ${CMAKE_SOURCE_DIR}/include)
TARGET_LINK_LIBRARIES(spacedSeeds ${CMAKE_SOURCE_DIR}/lib/libnthash.a unordered_dense::unordered_dense)
