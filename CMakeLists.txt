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

# compile flags
set(move_r_debug_flags "-Wall -Wextra -Wpedantic -ggdb")
set(move_r_release_flags "-w -march=native -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG           "${move_r_debug_flags} -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "${move_r_release_flags} -Ofast")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "${move_r_release_flags} -Os")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${move_r_debug_flags} -Ofast")

# options
option(MOVE_R_BUILD_CLI "Build cli programs" ON)
option(MOVE_R_BUILD_LZENDSA_CLI "Builds LZend-SA binaries" ON)
option(MOVE_R_BUILD_RLZSA_CLI "Builds RLZ-SA binaries" ON)
option(MOVE_R_BUILD_BENCH "Build the Move-r benchmark cli program" ON)
option(MOVE_R_BUILD_EXAMPLES "Build the example programs" ON)
option(MOVE_R_BUILD_TESTS "Build tests" ON)
option(MOVE_R_BUILD_INTERNAL_BENCH "Build programs for benchmarking Move-r's internal data structures" ON)
option(MOVE_R_USE_MALLOC_COUNT "Track memory usage via malloc_count (overrides malloc/free)" ON)

############################# move-r dependencies #############################

# openmp
find_package(OpenMP REQUIRED)

# ips4o
add_subdirectory("${CMAKE_SOURCE_DIR}/external/ips4o")
set_target_properties(ips4o_example PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

# ips2ra
set(IPS2RA_DISABLE_PARALLEL ON CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/ips2ra")

# libsais
set(LIBSAIS_USE_OPENMP ON CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/libsais/")

# malloc_count
add_library(malloc_count INTERFACE)
target_include_directories(malloc_count INTERFACE "${CMAKE_SOURCE_DIR}/external/malloc_count/")
if(MOVE_R_USE_MALLOC_COUNT)
  target_sources(malloc_count INTERFACE "${CMAKE_SOURCE_DIR}/external/malloc_count/malloc_count.c")
  target_link_libraries(malloc_count INTERFACE dl)
  target_compile_definitions(malloc_count INTERFACE MOVE_R_USE_MALLOC_COUNT=1)
else()
  target_sources(malloc_count INTERFACE "${CMAKE_SOURCE_DIR}/cmake/malloc_count_stub.c")
  target_compile_definitions(malloc_count INTERFACE MOVE_R_USE_MALLOC_COUNT=0)
endif()

# sdsl
set(SDSL_ROOT "${CMAKE_SOURCE_DIR}/external/sdsl-lite")
set(SDSL_INCLUDE_DIR "${SDSL_ROOT}/include")

# libdivsufsort (bundled with sdsl-lite).
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "" FORCE)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
add_subdirectory("${SDSL_ROOT}/external/libdivsufsort" "${CMAKE_BINARY_DIR}/external/libdivsufsort")
unset(CMAKE_POLICY_VERSION_MINIMUM)
set(DIVSUFSORT_INCLUDE_DIR "${CMAKE_BINARY_DIR}/external/libdivsufsort/include")
set_target_properties(divsufsort divsufsort64 PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(divsufsort INTERFACE "${DIVSUFSORT_INCLUDE_DIR}")
target_include_directories(divsufsort64 INTERFACE "${DIVSUFSORT_INCLUDE_DIR}")

# sdsl static library
file(GLOB sdsl_sources "${SDSL_ROOT}/lib/*.cpp")
add_library(sdsl STATIC ${sdsl_sources})
set_target_properties(sdsl PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(sdsl PRIVATE -w -DNDEBUG -fpermissive)
target_include_directories(sdsl PUBLIC "${SDSL_INCLUDE_DIR}")
target_include_directories(sdsl PRIVATE "${DIVSUFSORT_INCLUDE_DIR}")
target_link_libraries(sdsl PUBLIC divsufsort divsufsort64)

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
  OpenMP::OpenMP_CXX sdsl libsais ips4o ips2ra gtl sparse_map malloc_count
  emhash sux ordered rmq unordered_dense divsufsort divsufsort64
)

function(move_r_add_tool name src outdir)
  add_executable(${name} "${src}")
  target_link_libraries(${name} PRIVATE move_r ${ARGN})
  set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${outdir}")
endfunction()

############################# move-r cli #############################

if(MOVE_R_BUILD_CLI)
  foreach(tool move-r-build move-r-revert move-r-count move-r-locate move-r-query)
    move_r_add_tool(${tool} cli/move-r/${tool}.cpp "${CMAKE_BINARY_DIR}/cli/")
  endforeach()

  foreach(tool move-rb-build move-rb-revert move-rb-count move-rb-locate move-rb-query)
    move_r_add_tool(${tool} cli/move-rb/${tool}.cpp "${CMAKE_BINARY_DIR}/cli/")
  endforeach()

  move_r_add_tool(gen-patterns cli/tools/gen-patterns.cpp "${CMAKE_BINARY_DIR}/cli/")
  move_r_add_tool(print-data-structures cli/tools/print-data-structures.cpp "${CMAKE_BINARY_DIR}/cli/")
endif()

if(MOVE_R_BUILD_LZENDSA_CLI)
  foreach(tool lzendsa-build lzendsa-convert lzendsa-count lzendsa-locate lzendsa-random-access
               r-index-lzendsa-build r-index-lzendsa-count r-index-lzendsa-locate)
    move_r_add_tool(${tool} cli/lzendsa/${tool}.cpp "${CMAKE_BINARY_DIR}/cli/")
  endforeach()
endif()

if(MOVE_R_BUILD_RLZSA_CLI)
  foreach(tool rlzsa-build rlzsa-convert rlzsa-count rlzsa-locate rlzsa-random-access
               r-index-rlzsa-build r-index-rlzsa-count r-index-rlzsa-locate)
    move_r_add_tool(${tool} cli/rlzsa/${tool}.cpp "${CMAKE_BINARY_DIR}/cli/")
  endforeach()
endif()

############################# move-r examples #############################

if(MOVE_R_BUILD_EXAMPLES)
  foreach(ex move_data_structure move_r-queries move_r-retrieval move_r-store-load move-rb-queries)
    move_r_add_tool(example-${ex} examples/${ex}.cpp "${CMAKE_BINARY_DIR}/examples/")
  endforeach()
endif()

############################# move-r benchmark cli #############################

if(MOVE_R_BUILD_BENCH)
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

  move_r_add_tool(move-r-bench bench/move-r-bench/move-r-bench.cpp
    "${CMAKE_BINARY_DIR}/bench/" rcomp r-index-f r-index OnlineRlbwt)
endif()

############################# move_rb approximate-pattern-matching benchmark #############################

if(MOVE_R_BUILD_BENCH)
  set(columba_src "${CMAKE_SOURCE_DIR}/external/columba/src")
  set(columba_common_sources
    "${columba_src}/indexhelpers.cpp"
    "${columba_src}/bitparallelmatrix.cpp"
    "${columba_src}/searchstrategy.cpp"
    "${columba_src}/nucleotide.cpp"
    "${columba_src}/search.cpp"
    "${columba_src}/logger.cpp"
    "${columba_src}/parameters/alignparameters.cpp"
    "${columba_src}/indexinterface.cpp")

  function(move_r_add_columba_plugin name)
    add_library(${name} SHARED
      "${CMAKE_SOURCE_DIR}/bench/move-rb-bench/indexes/columba_capi.cpp" ${columba_common_sources} ${ARGN})
    target_include_directories(${name} PRIVATE
      "${columba_src}" "${CMAKE_SOURCE_DIR}/external/fmt/include" "${CMAKE_SOURCE_DIR}/external/parallel-hashmap"
      "${CMAKE_SOURCE_DIR}/include"
      "$<TARGET_PROPERTY:sdsl,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:gtl,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:sparse_map,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:ips4o,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:ips2ra,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:unordered_dense,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:ordered,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:sux,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:emhash,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:rmq,INTERFACE_INCLUDE_DIRECTORIES>"
      "$<TARGET_PROPERTY:malloc_count,INTERFACE_INCLUDE_DIRECTORIES>")
    target_compile_options(${name} PRIVATE -O3 -fno-fast-math -fno-strict-aliasing)
    target_link_libraries(${name} PRIVATE sdsl)
    set_target_properties(${name} PROPERTIES
      POSITION_INDEPENDENT_CODE ON LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bench/")
  endfunction()
  move_r_add_columba_plugin(columba_rlc_plugin
    "${columba_src}/bmove/bmove.cpp" "${columba_src}/bmove/moverepr.cpp")
  target_compile_definitions(columba_rlc_plugin PRIVATE ALPHABET=5 RUN_LENGTH_COMPRESSION PHI_MOVE FMT_HEADER_ONLY)
  move_r_add_columba_plugin(columba_fm_plugin
    "${columba_src}/fmindex/fmindex.cpp")
  target_compile_definitions(columba_fm_plugin PRIVATE ALPHABET=5 THIRTY_TWO FMT_HEADER_ONLY)

  add_library(br_index INTERFACE)
  target_include_directories(br_index INTERFACE "${CMAKE_SOURCE_DIR}/external/br-index/src")
  target_link_libraries(br_index INTERFACE sdsl divsufsort divsufsort64)

  move_r_add_tool(move-rb-bench bench/move-rb-bench/move-rb-bench.cpp
    "${CMAKE_BINARY_DIR}/bench/" br_index ${CMAKE_DL_LIBS})
  target_link_options(move-rb-bench PRIVATE -rdynamic)
  add_dependencies(move-rb-bench columba_rlc_plugin columba_fm_plugin)
  move_r_add_tool(move-rb-gen-queries bench/move-rb-bench/move-rb-gen-queries.cpp
    "${CMAKE_BINARY_DIR}/bench/")
endif()

############################# move-r internal bench #############################

if(MOVE_R_BUILD_INTERNAL_BENCH)
  move_r_add_tool(bench-int-rank-select bench/bench-int-rank-select.cpp
    "${CMAKE_BINARY_DIR}/bench/")
endif()

############################# move-r tests #############################

if(MOVE_R_BUILD_TESTS)
  add_subdirectory("${CMAKE_SOURCE_DIR}/external/googletest/")
  include(GoogleTest)
  enable_testing()

  set_target_properties(gtest gtest_main gmock gmock_main
    PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "")

  foreach(s move_r move_rb move_data_structure data_structures sa_index rlzsa_opt edit_distance_matrix)
    move_r_add_tool(test-${s} "${CMAKE_SOURCE_DIR}/tests/test_${s}.cpp"
      "${CMAKE_BINARY_DIR}/tests/" gtest_main)
  endforeach()

  gtest_discover_tests(test-move_r)
  gtest_discover_tests(test-move_rb)
  gtest_discover_tests(test-move_data_structure)
  gtest_discover_tests(test-data_structures)
  gtest_discover_tests(test-sa_index)
  gtest_discover_tests(test-rlzsa_opt)
  gtest_discover_tests(test-edit_distance_matrix)
endif()