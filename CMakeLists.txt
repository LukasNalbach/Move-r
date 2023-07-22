cmake_minimum_required(VERSION "3.16.3")
set (CMAKE_CXX_STANDARD 20)

find_package(Git REQUIRED)

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/libsais/LICENSE")
  execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
  execute_process(COMMAND cp -rf "${CMAKE_SOURCE_DIR}/fixed-files/*" "${CMAKE_SOURCE_DIR}")
endif()

project (move-r LANGUAGES C CXX)

cmake_policy(SET CMP0077 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
cmake_policy(SET CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(move_r_debug_flags "-Wall -Wextra -Wpedantic -ggdb") #-DDEBUG
set(move_r_release_flags "-w -march=native -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG           "${move_r_debug_flags} -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "${move_r_release_flags} -Ofast")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "${move_r_release_flags} -Os")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${move_r_debug_flags} -Ofast")

############################# move-r dependencies #############################

# openmp
find_package(OpenMP REQUIRED)

# ips4o
add_subdirectory(external/ips4o)
set_target_properties(ips4o_example PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

# concurrentqueue
add_subdirectory(external/concurrentqueue/)
include_directories(external/concurrentqueue/)

# libsais
set(LIBSAIS_USE_OPENMP ON)
add_subdirectory(external/libsais/)

# malloc_count
add_library(malloc_count INTERFACE)
target_sources(malloc_count INTERFACE external/malloc_count/malloc_count.c)
target_include_directories(malloc_count INTERFACE external/malloc_count/)
target_link_libraries(malloc_count INTERFACE dl)

# sdsl
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
add_subdirectory(external/sdsl-lite/)
set_target_properties(
  uninstall-sdsl uninstall gtest gtest_main divsufsort divsufsort64
  PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1
)
target_include_directories(sdsl PUBLIC ${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/include/)

# abseil
set(ABSL_PROPAGATE_CXX_STD ON CACHE BOOL "" FORCE)
add_subdirectory(external/abseil-cpp/)

############################# move-r #############################

# move-r
add_library(move_r INTERFACE)
target_include_directories(move_r INTERFACE include/)
target_link_libraries(move_r INTERFACE OpenMP::OpenMP_CXX sdsl libsais ips4o absl::btree malloc_count)

option(MOVE_R_BUILD_CLI "Build the cli programs" ON)
option(MOVE_R_BUILD_BENCH_CLI "Builds the benchmark cli program" OFF)
option(MOVE_R_BUILD_EXAMPLES "Build the example programs" OFF)

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
  add_subdirectory(external/Big-BWT/)

  ############################# other indexes #############################

  # rcomp
  add_library(rcomp INTERFACE)
  target_include_directories(rcomp INTERFACE external/rcomp/include/)

  # r-index-f
  add_library(r-index-f INTERFACE)
  target_include_directories(r-index-f INTERFACE
    external/r-index-f/include/common/
    external/r-index-f/include/r_index_f/
  )
  target_link_libraries(r-index-f INTERFACE sdsl)
  add_executable(pfp-thresholds external/pfp-thresholds/test/src/pfp_thresholds.cpp)
  target_sources(pfp-thresholds PRIVATE external/gsa-is/gsacak.c)
  target_include_directories(pfp-thresholds PRIVATE
    external/pfp-thresholds/include/common/
    external/pfp-thresholds/include/pfp/
    external/gsa-is/
  )
  target_link_libraries(pfp-thresholds PRIVATE sdsl malloc_count)
  set_target_properties(pfp-thresholds PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/external/pfp-thresholds/")

  # r-index-prezza
  add_library(r-index-prezza INTERFACE)
  target_include_directories(r-index-prezza INTERFACE external/r-index-prezza/internal/)
  target_link_libraries(r-index-prezza INTERFACE sdsl)

  # r-index-mun
  add_library(r-index-mun INTERFACE)
  target_include_directories(r-index-mun INTERFACE
    external/r-index-mun/internal/
    external/r-index-mun/
  )
  target_link_libraries(r-index-mun INTERFACE sdsl divsufsort divsufsort64)

  # OnlineRlbwt
  add_library(OnlineRlbwt INTERFACE)
  target_include_directories(OnlineRlbwt INTERFACE
    external/OnlineRlbwt/
    external/OnlineRlbwt/module/Basics/
    external/OnlineRlbwt/module/BTree/
  )

  # DYNAMIC
  add_library(DYNAMIC INTERFACE)
  target_include_directories(DYNAMIC INTERFACE
    external/DYNAMIC/include/
    external/hopscotch-map/include
  )

  ############################# move-r #############################

  add_executable(move-r-bench cli/move-r-bench.cpp)
  target_link_libraries(move-r-bench PRIVATE move_r rcomp r-index-f r-index-prezza r-index-mun OnlineRlbwt DYNAMIC)
  set_target_properties(move-r-bench PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/")
endif()