cmake_minimum_required(VERSION 3.10)

project(LogUP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)
include_directories(goldilocks/src)

# build goldilocks.a
file(GLOB DEP_SOURCES
    goldilocks/src/*.cpp
)
add_library(goldilocks STATIC ${DEP_SOURCES})
target_compile_options(goldilocks PRIVATE -O3 -mavx2 -MMD -MP -pthread -fopenmp)
find_package(OpenMP REQUIRED)
target_link_libraries(goldilocks PRIVATE
    OpenMP::OpenMP_CXX
    pthread
)

file(GLOB SOURCES
    src/*.cpp
)

# target test
add_executable(test test.cpp)
target_sources(test PRIVATE ${SOURCES})
target_compile_options(test PRIVATE -Wall -O3 -mavx2 -fopenmp)
target_link_libraries(test PRIVATE pthread gmp ssl crypto goldilocks)

# target bench
add_executable(bench bench.cpp)
target_sources(bench PRIVATE ${SOURCES})
target_compile_options(bench PRIVATE -Wall -O3 -mavx2 -fopenmp)
target_link_libraries(bench PRIVATE pthread gmp ssl crypto goldilocks)
