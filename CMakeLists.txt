cmake_minimum_required(VERSION "3.10.2")
set (CMAKE_CXX_STANDARD 17)

find_package(Git REQUIRED)

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/libsais/LICENSE")
  execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()

execute_process(COMMAND cp -rf "${CMAKE_SOURCE_DIR}/patched-files/*" "${CMAKE_SOURCE_DIR}")

project (move-r LANGUAGES C CXX)

cmake_policy(SET CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("Building in ${CMAKE_BUILD_TYPE} mode")

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
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/sdsl-lite/")
set_target_properties(
  uninstall-sdsl uninstall divsufsort divsufsort64
  PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1
)
target_include_directories(sdsl PUBLIC "${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/include/")

# abseil
set(ABSL_PROPAGATE_CXX_STD ON CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/abseil-cpp/")

############################# move-r #############################

# move-r
add_library(move_r INTERFACE)
target_include_directories(move_r INTERFACE "${CMAKE_SOURCE_DIR}/include/")
target_link_libraries(move_r INTERFACE OpenMP::OpenMP_CXX sdsl libsais ips4o absl::btree malloc_count concurrentqueue)

option(MOVE_R_BUILD_CLI "Build the cli programs" ON)
option(MOVE_R_BUILD_BENCH_CLI "Build the benchmark cli program" ON)
option(MOVE_R_BUILD_EXAMPLES "Build the example programs" OFF)
option(MOVE_R_BUILD_TESTS "Build tests" OFF)

############################# move-r cli #############################

if(MOVE_R_BUILD_CLI)
  add_executable(move-r-build cli/move-r-build.cpp)
  add_executable(move-r-revert cli/move-r-revert.cpp)
  add_executable(move-r-count cli/move-r-count.cpp)
  add_executable(move-r-locate cli/move-r-locate.cpp)
  add_executable(move-r-patterns cli/move-r-patterns.cpp)

  target_link_libraries(move-r-build PRIVATE move_r)
  target_link_libraries(move-r-revert PRIVATE move_r)
  target_link_libraries(move-r-count PRIVATE move_r)
  target_link_libraries(move-r-locate PRIVATE move_r)
  target_link_libraries(move-r-patterns PRIVATE move_r)
  
  set_target_properties(
    move-r-build move-r-revert move-r-count
    move-r-locate move-r-patterns
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/"
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
  # Big-BWT
  add_subdirectory("${CMAKE_SOURCE_DIR}/external/Big-BWT/")

  ############################# other indexes #############################

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
  target_include_directories(r-index INTERFACE "${CMAKE_SOURCE_DIR}/external/r-index/internal/")
  target_link_libraries(r-index INTERFACE sdsl)

  # r-index-bigbwt
  add_library(r-index-bigbwt INTERFACE)
  target_include_directories(r-index-bigbwt INTERFACE
    "${CMAKE_SOURCE_DIR}/external/r-index-mun/internal/"
    "${CMAKE_SOURCE_DIR}/external/r-index-mun/"
  )
  target_link_libraries(r-index-bigbwt INTERFACE sdsl divsufsort divsufsort64)

  # OnlineRlbwt
  add_library(OnlineRlbwt INTERFACE)
  target_include_directories(OnlineRlbwt INTERFACE
    "${CMAKE_SOURCE_DIR}/external/OnlineRlbwt/"
    "${CMAKE_SOURCE_DIR}/external/OnlineRlbwt/module/Basics/"
    "${CMAKE_SOURCE_DIR}/external/OnlineRlbwt/module/BTree/"
  )

  # DYNAMIC
  add_library(DYNAMIC INTERFACE)
  target_include_directories(DYNAMIC INTERFACE
    "${CMAKE_SOURCE_DIR}/external/DYNAMIC/include/"
    "${CMAKE_SOURCE_DIR}/external/hopscotch-map/include"
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

  add_executable(move-r-bench cli/move-r-bench/move-r-bench.cpp)
  target_link_libraries(move-r-bench PRIVATE move_r rcomp r-index-f r-index r-index-bigbwt OnlineRlbwt DYNAMIC)
  set_target_properties(move-r-bench PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/")
endif()

############################# move-r tests #############################

if(MOVE_R_BUILD_TESTS)
  enable_testing()
  include(GoogleTest)
  include(Dart)

  add_executable(test-move_r "${CMAKE_SOURCE_DIR}/tests/test_move_r.cpp")
  add_executable(test-move_data_structure "${CMAKE_SOURCE_DIR}/tests/test_move_data_structure.cpp")

  target_link_libraries(test-move_r PRIVATE gtest_main move_r)
  target_link_libraries(test-move_data_structure PRIVATE gtest_main move_r)

  set_target_properties(
    test-move_r test-move_data_structure
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/tests/"
  )

  gtest_discover_tests(test-move_r)
  gtest_discover_tests(test-move_data_structure)
endif()