cmake_minimum_required (VERSION 3.4)
project (gPartitioner C CXX)
set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC -lgomp -lrt -fopenmp -W -Wall -pedantic -O3 -Ofast -march=native")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3")


# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

add_subdirectory(lib/kProcessor)
set(kProcessor_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/lib/kProcessor/include")

# Add include directories.
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/MQF/include")
include_directories (${kProcessor_INCLUDE_PATH})


# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------


#include_directories("${PROJECT_SOURCE_DIR}/includes")


add_executable (genesIndexing genes_indexing.cpp)
target_link_libraries (genesIndexing kProcessor pthread z)
target_include_directories(genesIndexing INTERFACE ${kProcessor_INCLUDE_PATH})
