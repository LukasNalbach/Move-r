cmake_minimum_required(VERSION "3.10.2")
set(CMAKE_CXX_STANDARD 20)

find_package(Git REQUIRED)

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/libsais/LICENSE")
  message(STATUS "Downloading the git submodules")
  execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

  message(STATUS "Applying the patched files to the submodules")
  file(COPY "${CMAKE_SOURCE_DIR}/patched-files/external" DESTINATION "${CMAKE_SOURCE_DIR}")
endif()

project (move-r LANGUAGES C CXX)

cmake_policy(SET CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("building in ${CMAKE_BUILD_TYPE} mode")

# options
option(MOVE_R_BUILD_MOVE_R_CLI  "Build the move-r binaries"  ON)
option(MOVE_R_BUILD_MOVE_RB_CLI "Build the move-rb cli tools " ON)
option(MOVE_R_BUILD_LZENDSA_CLI "Build LZend-SA binaries" ON)
option(MOVE_R_BUILD_RLZSA_CLI "Build RLZ-SA binaries" ON)
option(MOVE_R_BUILD_MOVE_R_BENCH  "Build move-r-bench" ON)
option(MOVE_R_BUILD_MOVE_RB_BENCH "Build move-rb-bench-apm/-ext & move-rb-gen-apm-/-ext-queries" ON)
option(MOVE_R_BUILD_EXAMPLES "Build the example programs" ON)
option(MOVE_R_BUILD_TESTS "Build tests" ON)
option(MOVE_R_BUILD_INTERNAL_BENCH "Build programs for benchmarking Move-r's internal data structures" ON)
option(MOVE_R_USE_MALLOC_COUNT "Track memory usage via malloc_count (overrides malloc/free)" ON)
option(MOVE_R_MARCH_NATIVE "Tune the build for the host CPU (-march=native); turn OFF for portable binaries" ON)
option(MOVE_R_ENABLE_LTO "Enable link-time optimization (LTO/IPO) for move-r targets" ON)

if(WIN32)
  set(MOVE_R_BUILD_MOVE_R_BENCH OFF CACHE BOOL "" FORCE)
  set(MOVE_R_BUILD_MOVE_RB_BENCH OFF CACHE BOOL "" FORCE)
endif()

# compile flags
if(MOVE_R_MARCH_NATIVE)
  set(move_r_arch_flags "-march=native")
else()
  set(move_r_arch_flags "")
endif()

add_compile_options(-pthread)

set(move_r_debug_flags "-Wall -Wextra -Wpedantic -ggdb")
set(move_r_release_flags "-w ${move_r_arch_flags} -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG           "${move_r_debug_flags} -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "${move_r_release_flags} -O3")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "${move_r_release_flags} -Os")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${move_r_debug_flags} -O3")

if(WIN32)
  include_directories(BEFORE "${CMAKE_SOURCE_DIR}/cmake/winshim")
  add_compile_definitions(_REENTRANT)
  set(IPS4O_OPTIMIZE_FOR_NATIVE OFF CACHE BOOL "" FORCE)
  set(IPS2RA_OPTIMIZE_FOR_NATIVE OFF CACHE BOOL "" FORCE)
  set(IPS4O_DISABLE_PARALLEL ON CACHE BOOL "" FORCE)

  if(NOT MINGW)
    add_compile_definitions(MSVC_COMPILER _HAS_DEPRECATED_RESULT_OF=1)
  else()
    get_filename_component(_mingw_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(_mingw_root "${_mingw_bin}" DIRECTORY)
    if(EXISTS "${_mingw_root}/lib/libgomp.a")
      set(OpenMP_gomp_LIBRARY "${_mingw_root}/lib/libgomp.a" CACHE FILEPATH "" FORCE)
    endif()
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-mcx16)
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    get_filename_component(_llvm_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(_llvm_root "${_llvm_bin}" DIRECTORY)
    if(EXISTS "${_llvm_root}/lib/libomp.lib")
      foreach(_lang C CXX)
        set(OpenMP_${_lang}_FLAGS "-fopenmp=libomp" CACHE STRING "" FORCE)
        set(OpenMP_${_lang}_LIB_NAMES "libomp" CACHE STRING "" FORCE)
      endforeach()
      set(OpenMP_libomp_LIBRARY "${_llvm_root}/lib/libomp.lib" CACHE FILEPATH "" FORCE)
      if(EXISTS "${_llvm_root}/bin/libomp.dll")
        list(APPEND MOVE_R_WIN_RUNTIME_DLLS "${_llvm_root}/bin/libomp.dll")
      endif()
    endif()
  endif()
  
  set(MOVE_R_WINDOWS_TBB_ROOT "" CACHE PATH "Root of an extracted oneAPI TBB (Windows) release")
  if(MOVE_R_WINDOWS_TBB_ROOT)
    include_directories(BEFORE "${MOVE_R_WINDOWS_TBB_ROOT}/include")
    find_library(MOVE_R_WIN_TBB_LIB NAMES tbb12 tbb
      PATHS "${MOVE_R_WINDOWS_TBB_ROOT}/lib/intel64/vc14" "${MOVE_R_WINDOWS_TBB_ROOT}/lib" NO_DEFAULT_PATH)
    file(GLOB _tbb_dll "${MOVE_R_WINDOWS_TBB_ROOT}/redist/intel64/vc14/tbb12.dll"
                       "${MOVE_R_WINDOWS_TBB_ROOT}/redist/intel64/vc14/tbb.dll"
                       "${MOVE_R_WINDOWS_TBB_ROOT}/bin/*tbb12*.dll"
                       "${MOVE_R_WINDOWS_TBB_ROOT}/bin/*tbb.dll")
    if(_tbb_dll)
      list(APPEND MOVE_R_WIN_RUNTIME_DLLS ${_tbb_dll})
    endif()
  else()
    message(WARNING "MOVE_R_WINDOWS_TBB_ROOT is not set: ips4o's parallel sort needs Intel "
      "oneTBB on Windows. Download an 'oneapi-tbb-*-win' release and set MOVE_R_WINDOWS_TBB_ROOT.")
  endif()
endif()

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

set(BIGBWT_DIR "${CMAKE_SOURCE_DIR}/external/Big-BWT")
if(NOT WIN32 AND (EXISTS "${BIGBWT_DIR}/makefile" OR EXISTS "${BIGBWT_DIR}/Makefile"))
  add_custom_target(bigbwt ALL COMMAND make WORKING_DIRECTORY "${BIGBWT_DIR}"
    COMMENT "Building bundled Big-BWT (bigbwt)")
  add_compile_definitions(MOVE_R_BIGBWT_DIR="${BIGBWT_DIR}")
endif()

add_library(malloc_count INTERFACE)
target_include_directories(malloc_count INTERFACE "${CMAKE_SOURCE_DIR}/external/malloc_count/")
if(MOVE_R_USE_MALLOC_COUNT AND NOT WIN32)
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
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
add_subdirectory("${SDSL_ROOT}/external/libdivsufsort" "${CMAKE_BINARY_DIR}/external/libdivsufsort")
unset(CMAKE_POLICY_VERSION_MINIMUM)
set(CMAKE_WARN_DEPRECATED ON CACHE BOOL "" FORCE)
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
target_compile_options(move_r INTERFACE -fno-strict-aliasing)
target_link_libraries(move_r INTERFACE
  OpenMP::OpenMP_CXX sdsl libsais ips4o ips2ra gtl sparse_map malloc_count
  emhash sux ordered rmq unordered_dense divsufsort divsufsort64
)

if(WIN32 AND MOVE_R_WIN_TBB_LIB)
  target_link_libraries(move_r INTERFACE "${MOVE_R_WIN_TBB_LIB}")
endif()

if(MINGW)
  target_link_libraries(move_r INTERFACE atomic dl)
  target_link_options(move_r INTERFACE -static)
endif()

if(MOVE_R_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT MOVE_R_IPO_SUPPORTED OUTPUT MOVE_R_IPO_ERROR LANGUAGES CXX)
  if(NOT MOVE_R_IPO_SUPPORTED)
    message(WARNING "MOVE_R_ENABLE_LTO is ON but LTO/IPO is not supported by the toolchain: ${MOVE_R_IPO_ERROR}")
  endif()
endif()

function(move_r_add_tool name src outdir)
  add_executable(${name} "${src}")
  target_link_libraries(${name} PRIVATE move_r ${ARGN})
  set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${outdir}")
  if(MOVE_R_ENABLE_LTO AND MOVE_R_IPO_SUPPORTED)
    set_target_properties(${name} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
  endif()
  foreach(_dll ${MOVE_R_WIN_RUNTIME_DLLS})
    add_custom_command(TARGET ${name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_dll}" "${outdir}")
  endforeach()
endfunction()

############################# move-r cli #############################

if(MOVE_R_BUILD_MOVE_R_CLI)
  foreach(tool move-r-build move-r-revert move-r-count move-r-locate move-r-query)
    move_r_add_tool(${tool} cli/move-r/${tool}.cpp "${CMAKE_BINARY_DIR}/cli/")
  endforeach()
endif()

if(MOVE_R_BUILD_MOVE_RB_CLI)
  foreach(tool move-rb-build move-rb-revert move-rb-count move-rb-locate move-rb-query)
    move_r_add_tool(${tool} cli/move-rb/${tool}.cpp "${CMAKE_BINARY_DIR}/cli/")
  endforeach()
endif()

if(MOVE_R_BUILD_MOVE_R_CLI OR MOVE_R_BUILD_MOVE_RB_CLI)
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
  foreach(ex move-data-structure move-r-queries move-r-retrieval move-r-store-load move-rb-queries)
    move_r_add_tool(example-${ex} examples/${ex}.cpp "${CMAKE_BINARY_DIR}/examples/")
  endforeach()
endif()

############################# move-r benchmark cli #############################

if(MOVE_R_BUILD_MOVE_R_BENCH)
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
  set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
  add_subdirectory("${GRLBWT_ROOT}/external/cdt/")
  add_subdirectory("${GRLBWT_ROOT}/external/bioparsers/")
  set(CMAKE_WARN_DEPRECATED ON CACHE BOOL "" FORCE)
  add_library(grlbwt STATIC
    "${GRLBWT_ROOT}/lib/exact_algo/exact_ind_phase.cpp"
    "${GRLBWT_ROOT}/lib/exact_algo/exact_par_phase.cpp"
    "${GRLBWT_ROOT}/lib/opt_algo/opt_ind_phase.cpp"
    "${GRLBWT_ROOT}/lib/opt_algo/opt_par_phase.cpp"
    "${CMAKE_SOURCE_DIR}/external/malloc_count/malloc_count.c"
  )
  target_link_libraries(grlbwt PUBLIC z pthread cdt dl sdsl biopar)
  target_include_directories(grlbwt PUBLIC
    "${GRLBWT_ROOT}/include/"
    "${GRLBWT_ROOT}/include/exact_algo/"
    "${GRLBWT_ROOT}/include/opt_algo/"
    "${GRLBWT_ROOT}/external/"
    "${CMAKE_SOURCE_DIR}/external/malloc_count/"
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

if(MOVE_R_BUILD_MOVE_RB_BENCH)
  set(columba_src "${CMAKE_SOURCE_DIR}/external/columba/src")
  set(MOVE_R_COLUMBA_BITS "both" CACHE STRING "columba index widths to build for the move-rb-bench-* tools: both, 32 or 64")
  set_property(CACHE MOVE_R_COLUMBA_BITS PROPERTY STRINGS both 32 64)
  if(NOT MOVE_R_COLUMBA_BITS MATCHES "^(both|32|64)$")
    message(FATAL_ERROR "MOVE_R_COLUMBA_BITS must be one of: both, 32, 64 (got '${MOVE_R_COLUMBA_BITS}')")
  endif()
  set(columba_bits_32 FALSE)
  set(columba_bits_64 FALSE)
  if(MOVE_R_COLUMBA_BITS STREQUAL "both" OR MOVE_R_COLUMBA_BITS STREQUAL "32")
    set(columba_bits_32 TRUE)
  endif()
  if(MOVE_R_COLUMBA_BITS STREQUAL "both" OR MOVE_R_COLUMBA_BITS STREQUAL "64")
    set(columba_bits_64 TRUE)
  endif()
  set(columba_common_sources
    "${columba_src}/indexhelpers.cpp"
    "${columba_src}/bitparallelmatrix.cpp"
    "${columba_src}/searchstrategy.cpp"
    "${columba_src}/nucleotide.cpp"
    "${columba_src}/search.cpp"
    "${columba_src}/logger.cpp"
    "${columba_src}/parameters/alignparameters.cpp"
    "${columba_src}/indexinterface.cpp")

  set(columba_includes
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

  function(move_r_add_columba_plugin name)
    add_library(${name} SHARED
      "${CMAKE_SOURCE_DIR}/bench/move-rb-bench/indexes/columba_capi.cpp" ${columba_common_sources} ${ARGN})
    target_include_directories(${name} PRIVATE ${columba_includes})
    target_compile_options(${name} PRIVATE -O3 -fno-fast-math -fno-strict-aliasing)
    target_link_libraries(${name} PRIVATE sdsl)
    set_target_properties(${name} PROPERTIES
      POSITION_INDEPENDENT_CODE ON LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bench/")
  endfunction()
  if(columba_bits_32)
    move_r_add_columba_plugin(columba_rlc_plugin_32
      "${columba_src}/bmove/bmove.cpp" "${columba_src}/bmove/moverepr.cpp")
    target_compile_definitions(columba_rlc_plugin_32 PRIVATE ALPHABET=5 RUN_LENGTH_COMPRESSION PHI_MOVE FMT_HEADER_ONLY THIRTY_TWO)
    move_r_add_columba_plugin(columba_fm_plugin_32
      "${columba_src}/fmindex/fmindex.cpp")
    target_compile_definitions(columba_fm_plugin_32 PRIVATE ALPHABET=5 FMT_HEADER_ONLY THIRTY_TWO)
  endif()
  if(columba_bits_64)
    move_r_add_columba_plugin(columba_rlc_plugin_64
      "${columba_src}/bmove/bmove.cpp" "${columba_src}/bmove/moverepr.cpp")
    target_compile_definitions(columba_rlc_plugin_64 PRIVATE ALPHABET=5 RUN_LENGTH_COMPRESSION PHI_MOVE FMT_HEADER_ONLY)
    move_r_add_columba_plugin(columba_fm_plugin_64
      "${columba_src}/fmindex/fmindex.cpp")
    target_compile_definitions(columba_fm_plugin_64 PRIVATE ALPHABET=5 FMT_HEADER_ONLY)
  endif()

  add_library(br_index INTERFACE)
  target_include_directories(br_index INTERFACE "${CMAKE_SOURCE_DIR}/external/br-index/src")
  target_link_libraries(br_index INTERFACE sdsl divsufsort divsufsort64)

  ############################# index-build CLI tools for the compared indexes #############################

  # br-index's own index builder
  add_executable(bri-build "${CMAKE_SOURCE_DIR}/external/br-index/src/bri-build.cpp")
  target_link_libraries(bri-build PRIVATE br_index)
  set_target_properties(bri-build PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/")

  set(columba_build_sources
    "${columba_src}/buildindex.cpp" "${columba_src}/logger.cpp"
    "${columba_src}/parameters/buildparameters.cpp" "${columba_src}/seqfile.cpp")

  set(columba_bigbwt_dir "${CMAKE_SOURCE_DIR}/external/columba/external/Big-BWT")
  set(COLUMBA_BIGBWT_USABLE FALSE)
  if(EXISTS "${columba_bigbwt_dir}/utils.h")
    set(COLUMBA_BIGBWT_USABLE TRUE)
  endif()

  find_package(ZLIB REQUIRED)

  function(move_r_add_columba_build name rlc bits)
    add_executable(${name} ${columba_build_sources})
    target_include_directories(${name} PRIVATE ${columba_includes})
    target_compile_definitions(${name} PRIVATE ALPHABET=5 FMT_HEADER_ONLY HAVE_ZLIB)
    if(rlc)
      target_sources(${name} PRIVATE "${columba_src}/util.cpp" "${columba_src}/bmove/moverepr.cpp")
      target_compile_definitions(${name} PRIVATE RUN_LENGTH_COMPRESSION PHI_MOVE)
      if(COLUMBA_BIGBWT_USABLE)
        target_compile_definitions(${name} PRIVATE BIG_BWT_USABLE)
        target_include_directories(${name} PRIVATE "${columba_bigbwt_dir}")
      endif()
    endif()
    if(bits EQUAL 32)
      target_compile_definitions(${name} PRIVATE THIRTY_TWO)
      target_link_libraries(${name} PRIVATE libsais)
    else()
      target_link_libraries(${name} PRIVATE divsufsort64)
    endif()
    target_compile_options(${name} PRIVATE -O3 -fno-fast-math -fno-strict-aliasing)
    target_link_libraries(${name} PRIVATE sdsl ZLIB::ZLIB malloc_count)
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/columba-internal")
  endfunction()

  if(columba_bits_32)
    move_r_add_columba_build(columba-build-fm-32   FALSE 32)
    move_r_add_columba_build(columba-rlc-build-32  TRUE  32)
  endif()
  if(columba_bits_64)
    move_r_add_columba_build(columba-build-fm-64   FALSE 64)
    move_r_add_columba_build(columba-rlc-build-64  TRUE  64)
  endif()

  function(move_r_add_columba_build_wrapper name variant32 variant64)
    set(v32 "${variant32}")
    set(v64 "${variant64}")
    if(NOT v32)
      set(v32 "${v64}")
    endif()
    if(NOT v64)
      set(v64 "${v32}")
    endif()
    add_executable(${name} "${CMAKE_SOURCE_DIR}/cli/tools/columba-build-dispatch.cpp")
    target_compile_definitions(${name} PRIVATE
      COLUMBA_BUILD_32="columba-internal/${v32}" COLUMBA_BUILD_64="columba-internal/${v64}")
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/")
    if(variant32)
      add_dependencies(${name} ${variant32})
    endif()
    if(variant64)
      add_dependencies(${name} ${variant64})
    endif()
  endfunction()
  set(columba_fm_build_32 "")
  set(columba_fm_build_64 "")
  set(columba_rlc_build_32 "")
  set(columba_rlc_build_64 "")
  if(columba_bits_32)
    set(columba_fm_build_32 columba-build-fm-32)
    set(columba_rlc_build_32 columba-rlc-build-32)
  endif()
  if(columba_bits_64)
    set(columba_fm_build_64 columba-build-fm-64)
    set(columba_rlc_build_64 columba-rlc-build-64)
  endif()
  move_r_add_columba_build_wrapper(columba-build     "${columba_fm_build_32}"  "${columba_fm_build_64}")
  move_r_add_columba_build_wrapper(columba-rlc-build "${columba_rlc_build_32}" "${columba_rlc_build_64}")

  set(columba_plugins "")
  if(columba_bits_32)
    list(APPEND columba_plugins columba_rlc_plugin_32 columba_fm_plugin_32)
  endif()
  if(columba_bits_64)
    list(APPEND columba_plugins columba_rlc_plugin_64 columba_fm_plugin_64)
  endif()

  move_r_add_tool(move-rb-bench-apm bench/move-rb-bench/move-rb-bench-apm.cpp
    "${CMAKE_BINARY_DIR}/bench/" br_index ${CMAKE_DL_LIBS})
  target_link_options(move-rb-bench-apm PRIVATE -rdynamic)
  add_dependencies(move-rb-bench-apm ${columba_plugins})
  move_r_add_tool(move-rb-gen-apm-queries bench/move-rb-bench/move-rb-gen-apm-queries.cpp
    "${CMAKE_BINARY_DIR}/bench/")

  move_r_add_tool(move-rb-bench-ext bench/move-rb-bench/move-rb-bench-ext.cpp
    "${CMAKE_BINARY_DIR}/bench/" br_index ${CMAKE_DL_LIBS})
  target_link_options(move-rb-bench-ext PRIVATE -rdynamic)
  add_dependencies(move-rb-bench-ext ${columba_plugins})
  move_r_add_tool(move-rb-gen-ext-queries bench/move-rb-bench/move-rb-gen-ext-queries.cpp
    "${CMAKE_BINARY_DIR}/bench/")
endif()

############################# move-r tests #############################

if(MOVE_R_BUILD_TESTS)
  add_subdirectory("${CMAKE_SOURCE_DIR}/external/googletest/")
  include(GoogleTest)
  enable_testing()

  set_target_properties(gtest gtest_main gmock gmock_main
    PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "")

  foreach(s move-r move-rb move-data-structure data-structures csa-indexes rlzsa-opt edit-distance-matrix)
    move_r_add_tool("test-${s}" "${CMAKE_SOURCE_DIR}/tests/test-${s}.cpp"
      "${CMAKE_BINARY_DIR}/tests/" gtest_main)
    gtest_discover_tests("test-${s}")
  endforeach()
endif()