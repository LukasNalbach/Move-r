cmake_minimum_required(VERSION "3.10.2")
set(CMAKE_CXX_STANDARD 20)

find_package(Git REQUIRED)

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/libsais/LICENSE")
  execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()

project (move-r LANGUAGES C CXX)

cmake_policy(SET CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("building in ${CMAKE_BUILD_TYPE} mode")

set(move_r_debug_flags "-Wall -Wextra -Wpedantic -ggdb")
set(move_r_release_flags "-w -march=native -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG           "${move_r_debug_flags} -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "${move_r_release_flags} -Ofast")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "${move_r_release_flags} -Os")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${move_r_debug_flags} -Ofast")

############################# move-r dependencies #############################

# openmp
find_package(OpenMP REQUIRED)

# ips4o
add_subdirectory("${CMAKE_SOURCE_DIR}/external/ips4o")
set_target_properties(ips4o_example PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

# concurrentqueue
add_subdirectory("${CMAKE_SOURCE_DIR}/external/concurrentqueue/")

# libsais
set(LIBSAIS_USE_OPENMP ON CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/libsais/")

# malloc_count
add_library(malloc_count INTERFACE)
target_sources(malloc_count INTERFACE "${CMAKE_SOURCE_DIR}/external/malloc_count/malloc_count.c")
target_include_directories(malloc_count INTERFACE "${CMAKE_SOURCE_DIR}/external/malloc_count/")
target_link_libraries(malloc_count INTERFACE dl)

# sdsl
set(INCLUDE_DIRS "$ENV{HOME}/include" "$ENV{HOME}/local/include" "/usr/include" "/usr/local/include")
set(LIB_DIRS "$ENV{HOME}/lib" "$ENV{HOME}/local/lib /usr/lib" "/usr/local/lib")
find_path(SDSL_INCLUDE_DIR NAMES sdsl PATHS ${INCLUDE_DIRS})
find_library(SDSL_LIBRARY NAMES sdsl PATHS ${LIB_DIRS})
find_library(DIVSUFSORT_LIBRARY NAMES divsufsort PATHS ${LIB_DIRS})
find_library(DIVSUFSORT64_LIBRARY NAMES divsufsort64 PATHS ${LIB_DIRS})

if(NOT SDSL_INCLUDE_DIR OR NOT SDSL_LIBRARY)
  message(FATAL_ERROR "SDSL library not found")
endif()

add_library(sdsl STATIC IMPORTED GLOBAL)
add_library(divsufsort STATIC IMPORTED GLOBAL)
add_library(divsufsort64 STATIC IMPORTED GLOBAL)
set_target_properties(sdsl PROPERTIES IMPORTED_LOCATION "${SDSL_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES "${SDSL_INCLUDE_DIR}")
set_target_properties(divsufsort PROPERTIES IMPORTED_LOCATION "${DIVSUFSORT_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES "${SDSL_INCLUDE_DIR}")
set_target_properties(divsufsort64 PROPERTIES IMPORTED_LOCATION "${DIVSUFSORT64_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES "${SDSL_INCLUDE_DIR}")

# gtl
add_subdirectory("${CMAKE_SOURCE_DIR}/external/gtl/")

# sparse-map
add_subdirectory("${CMAKE_SOURCE_DIR}/external/sparse-map/")

# emhash
add_library(emhash INTERFACE)
target_include_directories(emhash INTERFACE "${CMAKE_SOURCE_DIR}/external/emhash/")

# sux
add_library(sux INTERFACE)
target_include_directories(sux INTERFACE "${CMAKE_SOURCE_DIR}/external/sux/")

# ordered
add_subdirectory("${CMAKE_SOURCE_DIR}/external/ordered/")

# rmq
add_library(rmq INTERFACE)
target_include_directories(rmq INTERFACE "${CMAKE_SOURCE_DIR}/external/rmq/include/")

# unordered_dense
add_subdirectory("${CMAKE_SOURCE_DIR}/external/unordered_dense/")

############################# move-r #############################

# move-r
add_library(move_r INTERFACE)
target_include_directories(move_r INTERFACE "${CMAKE_SOURCE_DIR}/include/")
target_link_libraries(move_r INTERFACE
  OpenMP::OpenMP_CXX sdsl libsais ips4o gtl sparse_map malloc_count
  concurrentqueue emhash sux ordered rmq unordered_dense divsufsort divsufsort64
)

option(MOVE_R_BUILD_CLI "Build cli programs" ON)
option(MOVE_R_BUILD_LZENDSA_CLI "Builds Rlzsa binaries" ON)
option(MOVE_R_BUILD_RLZSA_CLI "Builds Rlzsa binaries" ON)
option(MOVE_R_BUILD_BENCH_CLI "Build the benchmark cli program" ON)
option(MOVE_R_BUILD_EXAMPLES "Build the example programs" ON)
option(MOVE_R_BUILD_TESTS "Build tests" ON)
option(MOVE_R_BUILD_BENCH "Build program for benchmarking internal data structures" ON)

############################# move-r cli #############################

if(MOVE_R_BUILD_CLI)
  add_executable(move-r-build cli/move-r/move-r-build.cpp)
  add_executable(move-r-revert cli/move-r/move-r-revert.cpp)
  add_executable(move-r-count cli/move-r/move-r-count.cpp)
  add_executable(move-r-locate cli/move-r/move-r-locate.cpp)
  add_executable(move-r-patterns cli/move-r/move-r-patterns.cpp)

  add_executable(move-rb-build cli/move-rb/move-rb-build.cpp)
  add_executable(move-rb-count cli/move-rb/move-rb-count.cpp)
  add_executable(move-rb-locate cli/move-rb/move-rb-locate.cpp)
  add_executable(move-rb-patterns cli/move-rb/move-rb-patterns.cpp)

  target_link_libraries(move-r-build PRIVATE move_r)
  target_link_libraries(move-rb-build PRIVATE move_r)
  target_link_libraries(move-r-revert PRIVATE move_r)
  target_link_libraries(move-r-count PRIVATE move_r)
  target_link_libraries(move-r-locate PRIVATE move_r)
  target_link_libraries(move-r-patterns PRIVATE move_r)

  target_link_libraries(move-rb-build PRIVATE move_r)
  target_link_libraries(move-rb-count PRIVATE move_r)
  target_link_libraries(move-rb-locate PRIVATE move_r)
  target_link_libraries(move-rb-patterns PRIVATE move_r)

  set_target_properties(
    move-r-build move-r-revert move-r-count move-r-locate move-r-patterns
    move-rb-build move-rb-count move-rb-locate move-rb-patterns
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/"
  )
endif()

if(MOVE_R_BUILD_LZENDSA_CLI)
  add_executable(lzendsa-build cli/lzendsa/lzendsa-build.cpp)
  add_executable(lzendsa-convert cli/lzendsa/lzendsa-convert.cpp)
  add_executable(lzendsa-count cli/lzendsa/lzendsa-count.cpp)
  add_executable(lzendsa-locate cli/lzendsa/lzendsa-locate.cpp)
  add_executable(lzendsa-random-access cli/lzendsa/lzendsa-random-access.cpp)
  add_executable(r-index-lzendsa-build cli/lzendsa/r-index-lzendsa-build.cpp)
  add_executable(r-index-lzendsa-count cli/lzendsa/r-index-lzendsa-count.cpp)
  add_executable(r-index-lzendsa-locate cli/lzendsa/r-index-lzendsa-locate.cpp)

  target_link_libraries(lzendsa-build PRIVATE move_r)
  target_link_libraries(lzendsa-convert PRIVATE move_r)
  target_link_libraries(lzendsa-count PRIVATE move_r)
  target_link_libraries(lzendsa-locate PRIVATE move_r)
  target_link_libraries(lzendsa-random-access PRIVATE move_r)
  target_link_libraries(r-index-lzendsa-build PRIVATE move_r)
  target_link_libraries(r-index-lzendsa-count PRIVATE move_r)
  target_link_libraries(r-index-lzendsa-locate PRIVATE move_r)

  set_target_properties(
    lzendsa-build lzendsa-convert lzendsa-count lzendsa-locate lzendsa-random-access
    r-index-lzendsa-build r-index-lzendsa-count r-index-lzendsa-locate
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/"
  )
endif()

if(MOVE_R_BUILD_RLZSA_CLI)
  add_executable(rlzsa-build cli/rlzsa/rlzsa-build.cpp)
  add_executable(rlzsa-convert cli/rlzsa/rlzsa-convert.cpp)
  add_executable(rlzsa-count cli/rlzsa/rlzsa-count.cpp)
  add_executable(rlzsa-locate cli/rlzsa/rlzsa-locate.cpp)
  add_executable(rlzsa-random-access cli/rlzsa/rlzsa-random-access.cpp)
  add_executable(r-index-rlzsa-build cli/rlzsa/r-index-rlzsa-build.cpp)
  add_executable(r-index-rlzsa-count cli/rlzsa/r-index-rlzsa-count.cpp)
  add_executable(r-index-rlzsa-locate cli/rlzsa/r-index-rlzsa-locate.cpp)

  target_link_libraries(rlzsa-build PRIVATE move_r)
  target_link_libraries(rlzsa-convert PRIVATE move_r)
  target_link_libraries(rlzsa-count PRIVATE move_r)
  target_link_libraries(rlzsa-locate PRIVATE move_r)
  target_link_libraries(rlzsa-random-access PRIVATE move_r)
  target_link_libraries(r-index-rlzsa-build PRIVATE move_r)
  target_link_libraries(r-index-rlzsa-count PRIVATE move_r)
  target_link_libraries(r-index-rlzsa-locate PRIVATE move_r)
  
  set_target_properties(
    rlzsa-build rlzsa-convert rlzsa-count rlzsa-locate rlzsa-random-access
    r-index-rlzsa-build r-index-rlzsa-count r-index-rlzsa-locate
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/"
  )
endif()

############################# move-r examples #############################

if(MOVE_R_BUILD_EXAMPLES)
  add_executable(example-move_data_structure examples/move_data_structure.cpp)
  add_executable(example-move_r-queries examples/move_r-queries.cpp)
  add_executable(example-move_r-retrieval examples/move_r-retrieval.cpp)
  add_executable(example-move_r-store-load examples/move_r-store-load.cpp)

  target_link_libraries(example-move_data_structure PRIVATE move_r)
  target_link_libraries(example-move_r-queries PRIVATE move_r)
  target_link_libraries(example-move_r-retrieval PRIVATE move_r)
  target_link_libraries(example-move_r-store-load PRIVATE move_r)
  
  set_target_properties(
    example-move_data_structure example-move_r-queries
    example-move_r-retrieval example-move_r-store-load
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples/"
  )
endif()

############################# move-r benchmark cli #############################

if(MOVE_R_BUILD_BENCH_CLI)
  # rcomp
  add_library(rcomp INTERFACE)
  target_include_directories(rcomp INTERFACE "${CMAKE_SOURCE_DIR}/external/rcomp/include/")

  # r-index-f
  add_library(r-index-f INTERFACE)
  target_include_directories(r-index-f INTERFACE
    "${CMAKE_SOURCE_DIR}/external/r-index-f/include/common/"
    "${CMAKE_SOURCE_DIR}/external/r-index-f/include/r_index_f/"
  )
  target_link_libraries(r-index-f INTERFACE sdsl)
  add_executable(pfp-thresholds "${CMAKE_SOURCE_DIR}/external/pfp-thresholds/test/src/pfp_thresholds.cpp")
  target_sources(pfp-thresholds PRIVATE "${CMAKE_SOURCE_DIR}/external/gsa-is/gsacak.c")
  target_include_directories(pfp-thresholds PRIVATE
    "${CMAKE_SOURCE_DIR}/external/pfp-thresholds/include/common/"
    "${CMAKE_SOURCE_DIR}/external/pfp-thresholds/include/pfp/"
    "${CMAKE_SOURCE_DIR}/external/gsa-is/"
  )
  target_link_libraries(pfp-thresholds PRIVATE sdsl malloc_count)
  set_target_properties(pfp-thresholds PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/external/pfp-thresholds/")

  # r-index
  add_library(r-index INTERFACE)
  target_include_directories(r-index INTERFACE
    "${CMAKE_SOURCE_DIR}/external/r-index/internal/"
    "${CMAKE_SOURCE_DIR}/external/r-index/"
  )
  target_link_libraries(r-index INTERFACE sdsl divsufsort divsufsort64)

  # OnlineRlbwt
  add_library(OnlineRlbwt INTERFACE)
  target_include_directories(OnlineRlbwt INTERFACE
    "${CMAKE_SOURCE_DIR}/external/OnlineRlbwt/"
    "${CMAKE_SOURCE_DIR}/external/OnlineRlbwt/module/Basics/"
    "${CMAKE_SOURCE_DIR}/external/OnlineRlbwt/module/BTree/"
  )

  # grlBWT
  set(GRLBWT_ROOT "${CMAKE_SOURCE_DIR}/external/grlBWT/")
  add_subdirectory("${GRLBWT_ROOT}/external/cdt/")
  add_subdirectory("${GRLBWT_ROOT}/external/bioparsers/")
  add_library(grlbwt STATIC
    "${GRLBWT_ROOT}/lib/exact_algo/exact_ind_phase.cpp"
    "${GRLBWT_ROOT}/lib/exact_algo/exact_par_phase.cpp"
    "${GRLBWT_ROOT}/lib/opt_algo/opt_ind_phase.cpp"
    "${GRLBWT_ROOT}/lib/opt_algo/opt_par_phase.cpp"
    "${GRLBWT_ROOT}/external/malloc_count-master/malloc_count.c"
  )
  target_link_libraries(grlbwt PUBLIC z pthread cdt dl sdsl biopar)
  target_include_directories(grlbwt PUBLIC
    "${GRLBWT_ROOT}/include/"
    "${GRLBWT_ROOT}/include/exact_algo/"
    "${GRLBWT_ROOT}/include/opt_algo/"
    "${GRLBWT_ROOT}/external/"
    "${GRLBWT_ROOT}/external/malloc_count-master/"
  )
  add_executable(grlbwt-cli "${GRLBWT_ROOT}/main.cpp")
  add_executable(grlbwt2rle "${GRLBWT_ROOT}/scripts/grlbwt2rle.cpp")
  target_link_libraries(grlbwt-cli PRIVATE grlbwt)
  target_link_libraries(grlbwt2rle PRIVATE grlbwt)
  set_target_properties(grlbwt-cli grlbwt2rle
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/external/grlBWT/"
  )

  # block_RLBWT
  if(NOT EXISTS "${CMAKE_BINARY_DIR}/external/block_RLBWT/")
    execute_process(COMMAND mkdir "${CMAKE_BINARY_DIR}/external/block_RLBWT/")
  endif()
  execute_process(
    COMMAND make "${CMAKE_SOURCE_DIR}/external/block_RLBWT/make_alphabet_header"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/external/block_RLBWT/"
  )
  execute_process(COMMAND mv
    "${CMAKE_SOURCE_DIR}/external/block_RLBWT/make_alphabet_header"
    "${CMAKE_BINARY_DIR}/external/block_RLBWT/make_alphabet_header"
  )

  ############################# move-r #############################

  add_executable(move-r-bench cli/move-r/move-r-bench/move-r-bench.cpp)
  target_link_libraries(move-r-bench PRIVATE move_r rcomp r-index-f r-index OnlineRlbwt)
  set_target_properties(move-r-bench PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/")
endif()

############################# move-r internal bench #############################

if(MOVE_R_BUILD_BENCH)
  add_executable(bench-int-rank-select bench/bench_int_rank_select.cpp)
  target_link_libraries(bench-int-rank-select PRIVATE move_r)
  set_target_properties(bench-int-rank-select PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bench/")
endif()

############################# move-r tests #############################

if(MOVE_R_BUILD_TESTS)
  add_subdirectory("${CMAKE_SOURCE_DIR}/external/googletest/")
  include(GoogleTest)
  enable_testing()

  add_executable(test-move_r "${CMAKE_SOURCE_DIR}/tests/test_move_r.cpp")
  add_executable(test-move_r_int "${CMAKE_SOURCE_DIR}/tests/test_move_r_int.cpp")
  add_executable(test-move_rb "${CMAKE_SOURCE_DIR}/tests/test_move_rb.cpp")
  add_executable(test-move_data_structure "${CMAKE_SOURCE_DIR}/tests/test_move_data_structure.cpp")

  target_link_libraries(test-move_r PRIVATE gtest_main move_r)
  target_link_libraries(test-move_r_int PRIVATE gtest_main move_r)
  target_link_libraries(test-move_rb PRIVATE gtest_main move_r)
  target_link_libraries(test-move_data_structure PRIVATE gtest_main move_r)

  set_target_properties(
    test-move_r test-move_r_int test-move_rb test-move_data_structure
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/tests/"
  )

  gtest_discover_tests(test-move_r)
  gtest_discover_tests(test-move_r_int)
  gtest_discover_tests(test-move_data_structure)
endif()