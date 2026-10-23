cmake_minimum_required(VERSION 3.21)

# Move-r: the move-r and move-rb indexes, their command-line tools, benchmarks and tests. If Move-r is added by another project, only its library is built.
project (move-r LANGUAGES C CXX)

cmake_policy(SET CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

# string-utils (unless the project that added Move-r already provides it)
if(NOT TARGET string_utils)
  file(GLOB move_r_string_utils_contents "${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils/*")

  if(move_r_string_utils_contents STREQUAL "")
    find_package(Git REQUIRED)
    message(STATUS "Downloading the git submodule external/string-utils")
    execute_process(COMMAND "${GIT_EXECUTABLE}" submodule update --init -- external/string-utils
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  # the prefix of the options and environment variables (MOVE_R_MARCH_NATIVE, MOVE_R_TMPDIR, ...); they also apply
  # to the submodules added below
  if(PROJECT_IS_TOP_LEVEL)
    set(STRING_UTILS_PREFIX MOVE_R)
    include("${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils/cmake/settings.cmake")
  endif()

  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils" EXCLUDE_FROM_ALL)
endif()

# options: the command-line tools and the tests are ON if Move-r is the top-level project; the benchmarks are always
# OFF unless they are switched on, then their external indexes are downloaded and built too. The global build options
# are defined in external/string-utils/cmake/.
option(MOVE_R_BUILD_MOVE_R_CLI     "Build the move-r command line tools" ${PROJECT_IS_TOP_LEVEL})
option(MOVE_R_BUILD_MOVE_RB_CLI    "Build the move-rb command line tools" ${PROJECT_IS_TOP_LEVEL})
option(MOVE_R_BUILD_TOOLS_CLI      "Build the shared command line tools (gen-patterns, gen-repetitive, print-data-structures)" ${PROJECT_IS_TOP_LEVEL})
option(MOVE_R_BUILD_TESTS "Build tests" ${PROJECT_IS_TOP_LEVEL})

# benchmarks
option(MOVE_R_BUILD_MOVE_R_BENCH    "Build move-r-bench, move-r-bench-ops and the external indexes they measure against" OFF)
option(MOVE_R_BUILD_MOVE_RB_BENCH   "Build the move-rb benchmarks and the external indexes they measure against (columba, br-index)" OFF)

# download the missing submodules: the own components and sparse-map without their submodules (string-utils downloads
# the libraries it provides), and the external indexes of the benchmarks that are on with theirs
set(move_r_component_submodules external/move-data-structure external/rlz-ref external/pfp-rlbwt
  external/sparse-map)

set(move_r_third_party_submodules "")
if(MOVE_R_BUILD_MOVE_R_BENCH)
  list(APPEND move_r_third_party_submodules external/rcomp external/r-index-f external/pfp-thresholds
    external/gsa-is external/r-index external/OnlineRlbwt external/grlBWT external/block_RLBWT)
endif()
if(MOVE_R_BUILD_MOVE_RB_BENCH)
  list(APPEND move_r_third_party_submodules external/columba external/fmt external/parallel-hashmap
    external/br-index external/sassy)
endif()

string_utils_fetch_submodules("${CMAKE_CURRENT_SOURCE_DIR}" PATHS ${move_r_component_submodules})
string_utils_fetch_submodules("${CMAKE_CURRENT_SOURCE_DIR}" RECURSIVE PATHS ${move_r_third_party_submodules})

# the libraries of string-utils beyond its core that are used here (the external indexes of the benchmarks use all of
# sdsl)
string_utils_add_components(libsais gtl emhash)
if(MOVE_R_BUILD_MOVE_R_BENCH OR MOVE_R_BUILD_MOVE_RB_BENCH)
  string_utils_add_components(sdsl)
endif()

# apply the patches in patches/ to the downloaded submodules
string_utils_apply_patches("${CMAKE_CURRENT_SOURCE_DIR}")

# the benchmarks against external indexes (move-r-bench, the move-rb benchmarks) use POSIX (fork, dlopen, /proc) and are not built on Windows; the benchmarks of the operations of our indexes
# (move-r-bench-ops, move-rb-bench-ops) are
if(WIN32)
  set(move_r_external_benches OFF)
else()
  set(move_r_external_benches ON)
endif()

############################# move-r dependencies #############################

# OpenMP, sdsl, libsais, ips4o, ips2ra, malloc_count etc. are provided by string-utils

# sparse-map
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/sparse-map/" EXCLUDE_FROM_ALL)

# move-data-structure (unless the project that added Move-r already provides it)
if(NOT TARGET move_data_structure)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/move-data-structure" EXCLUDE_FROM_ALL)
endif()

# pfp-rlbwt (used by the construction mode _pfp)
if(NOT TARGET pfp_rlbwt)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/pfp-rlbwt" EXCLUDE_FROM_ALL)
endif()

# rlz-ref (unless the project that added Move-r already provides it)
if(NOT TARGET rlz_ref)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/rlz-ref" EXCLUDE_FROM_ALL)
endif()

############################# move-r #############################

# move-r
add_library(move_r INTERFACE)
target_include_directories(move_r INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_link_libraries(move_r INTERFACE string_utils libsais libsais40 gtl emhash move_data_structure rlz_ref
  sparse_map pfp_rlbwt)

# the prefix-free parsing of a text with an integer alphabet (see include/move_r/construction/pfp_integer.hpp)
add_library(move_r_pfp_integer STATIC src/pfp_integer.cpp)
target_include_directories(move_r_pfp_integer PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_link_libraries(move_r_pfp_integer PRIVATE string_utils sparse_map pfp_rlbwt)
target_link_libraries(move_r INTERFACE move_r_pfp_integer)

# explicit instantiations of all profiles of the indexes (the option MOVE_R_EXTERN_TEMPLATES, whose value string-utils
# stores in STRING_UTILS_EXTERN_TEMPLATES); each translation unit instantiates one group of members for all index
# types (see src/instantiations.cpp)
if(STRING_UTILS_EXTERN_TEMPLATES)
  # src/instantiations.cpp is compiled once per part
  set(move_r_instantiation_parts
    MOVE_R MOVE_R_CONSTRUCTION_TEXT MOVE_R_CONSTRUCTION_MOVE MOVE_R_CONSTRUCTION_SA_RLZ
    MOVE_RB_1 MOVE_RB_2 MOVE_RB_BUILD_1 MOVE_RB_BUILD_2 MOVE_RB_APM_1 MOVE_RB_APM_2
    SORT TEXT_OPS)
  set(move_r_instantiation_objects "")

  foreach(part IN LISTS move_r_instantiation_parts)
    string(TOLOWER "move_r_instantiations_${part}" part_target)
    add_library(${part_target} OBJECT src/instantiations.cpp)
    target_compile_definitions(${part_target} PRIVATE MOVE_R_INSTANTIATION_PART_${part})
    target_include_directories(${part_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/")
    target_link_libraries(${part_target} PRIVATE string_utils libsais libsais40 gtl emhash move_data_structure
      rlz_ref sparse_map pfp_rlbwt)
    list(APPEND move_r_instantiation_objects $<TARGET_OBJECTS:${part_target}>)
  endforeach()

  add_library(move_r_instantiations STATIC ${move_r_instantiation_objects})
  target_link_libraries(move_r_instantiations PRIVATE string_utils move_data_structure rlz_ref sparse_map
    move_r_pfp_integer pfp_rlbwt)
  target_link_libraries(move_r INTERFACE move_r_instantiations)
endif()

# adds an executable (see string_utils_add_executable) that is linked against move_r and the libraries in ARGN
function(move_r_add_tool name src outdir)
  string_utils_add_executable(${name} "${src}" "${outdir}" move_r ${ARGN})
endfunction()

############################# move-r cli #############################

# a command-line tool whose modes are implemented in translation units in cli/<dir>/ that are compiled in parallel
# (see cli/<dir>/modes.hpp)
function(move_r_add_cli name dir)
  list(TRANSFORM ARGN PREPEND "cli/${dir}/" OUTPUT_VARIABLE sources)
  list(TRANSFORM sources APPEND ".cpp")
  move_r_add_tool(${name} "${sources}" "${CMAKE_BINARY_DIR}/cli/")
endfunction()

if(MOVE_R_BUILD_MOVE_R_CLI)
  move_r_add_cli(move-r move-r main build query revert extract info)
endif()

if(MOVE_R_BUILD_MOVE_RB_CLI)
  move_r_add_cli(move-rb move-rb main build query revert extract info)
endif()

if(MOVE_R_BUILD_TOOLS_CLI)
  move_r_add_tool(gen-patterns cli/tools/gen-patterns.cpp "${CMAKE_BINARY_DIR}/cli/")
  move_r_add_tool(gen-repetitive cli/tools/gen-repetitive.cpp "${CMAKE_BINARY_DIR}/cli/")
  move_r_add_tool(fasta-to-text cli/tools/fasta-to-text.cpp "${CMAKE_BINARY_DIR}/cli/")
  move_r_add_tool(print-data-structures cli/tools/print-data-structures.cpp "${CMAKE_BINARY_DIR}/cli/")
  # compiled without optimization
  target_compile_options(print-data-structures PRIVATE -O0)
endif()

############################# internal benchmarks #############################

if(MOVE_R_BUILD_MOVE_R_BENCH)
  move_r_add_tool(move-r-bench-ops bench/move-r-bench-ops/move-r-bench-ops.cpp
    "${CMAKE_BINARY_DIR}/bench/")
endif()

if(MOVE_R_BUILD_MOVE_RB_BENCH)
  move_r_add_tool(move-rb-bench-ops bench/move-rb-bench-ops/move-rb-bench-ops.cpp
    "${CMAKE_BINARY_DIR}/bench/")
endif()

############################# move-r benchmark cli #############################

if(MOVE_R_BUILD_MOVE_R_BENCH AND move_r_external_benches)
  # rcomp
  add_library(rcomp INTERFACE)
  target_include_directories(rcomp INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/external/rcomp/include/")

  # r-index-f
  add_library(r-index-f INTERFACE)
  target_include_directories(r-index-f INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/external/r-index-f/include/common/"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/r-index-f/include/r_index_f/"
  )
  target_link_libraries(r-index-f INTERFACE sdsl)
  add_executable(pfp-thresholds "${CMAKE_CURRENT_SOURCE_DIR}/external/pfp-thresholds/test/src/pfp_thresholds.cpp")
  target_sources(pfp-thresholds PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/gsa-is/gsacak.c")
  target_include_directories(pfp-thresholds PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/external/pfp-thresholds/include/common/"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/pfp-thresholds/include/pfp/"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/gsa-is/"
  )
  target_link_libraries(pfp-thresholds PRIVATE sdsl malloc_count)
  set_target_properties(pfp-thresholds PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/external/pfp-thresholds/")

  # r-index
  add_library(r-index INTERFACE)
  target_include_directories(r-index INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/external/r-index/internal/"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/r-index/"
  )
  target_link_libraries(r-index INTERFACE sdsl divsufsort divsufsort64)

  # OnlineRlbwt
  add_library(OnlineRlbwt INTERFACE)
  target_include_directories(OnlineRlbwt INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/external/OnlineRlbwt/"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/OnlineRlbwt/module/Basics/"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/OnlineRlbwt/module/BTree/"
  )

  # grlBWT
  set(GRLBWT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/external/grlBWT/")
  set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
  add_subdirectory("${GRLBWT_ROOT}/external/cdt/")
  add_subdirectory("${GRLBWT_ROOT}/external/bioparsers/")
  set(CMAKE_WARN_DEPRECATED ON CACHE BOOL "" FORCE)
  add_library(grlbwt STATIC
    "${GRLBWT_ROOT}/lib/exact_algo/exact_ind_phase.cpp"
    "${GRLBWT_ROOT}/lib/exact_algo/exact_par_phase.cpp"
    "${GRLBWT_ROOT}/lib/opt_algo/opt_ind_phase.cpp"
    "${GRLBWT_ROOT}/lib/opt_algo/opt_par_phase.cpp"
    "${STRING_UTILS_EXTERNAL_DIR}/malloc_count/malloc_count.c"
  )
  target_link_libraries(grlbwt PUBLIC z pthread cdt dl sdsl biopar)
  target_include_directories(grlbwt PUBLIC
    "${GRLBWT_ROOT}/include/"
    "${GRLBWT_ROOT}/include/exact_algo/"
    "${GRLBWT_ROOT}/include/opt_algo/"
    "${GRLBWT_ROOT}/external/"
    "${STRING_UTILS_EXTERNAL_DIR}/malloc_count/"
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
    COMMAND make "${CMAKE_CURRENT_SOURCE_DIR}/external/block_RLBWT/make_alphabet_header"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/external/block_RLBWT/"
  )
  execute_process(COMMAND mv
    "${CMAKE_CURRENT_SOURCE_DIR}/external/block_RLBWT/make_alphabet_header"
    "${CMAKE_BINARY_DIR}/external/block_RLBWT/make_alphabet_header"
  )

  ############################# move-r #############################

  move_r_add_tool(move-r-bench bench/move-r-bench/move-r-bench.cpp
    "${CMAKE_BINARY_DIR}/bench/" rcomp r-index-f r-index OnlineRlbwt)
endif()

############################# move_rb approximate-pattern-matching benchmark #############################

if(MOVE_R_BUILD_MOVE_RB_BENCH AND move_r_external_benches)
  set(columba_src "${CMAKE_CURRENT_SOURCE_DIR}/external/columba/src")
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
    "${columba_src}" "${CMAKE_CURRENT_SOURCE_DIR}/external/fmt/include" "${CMAKE_CURRENT_SOURCE_DIR}/external/parallel-hashmap"
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "$<TARGET_PROPERTY:string_utils,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:move_data_structure,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:rlz_ref,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:pfp_rlbwt,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:sdsl,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:gtl,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:sparse_map,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:ips4o,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:ips2ra,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:sux,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:emhash,INTERFACE_INCLUDE_DIRECTORIES>"
    "$<TARGET_PROPERTY:malloc_count,INTERFACE_INCLUDE_DIRECTORIES>")

  function(move_r_add_columba_plugin name)
    add_library(${name} SHARED
      "${CMAKE_CURRENT_SOURCE_DIR}/bench/move-rb-bench/indexes/columba_capi.cpp" ${columba_common_sources} ${ARGN})
    target_include_directories(${name} SYSTEM PRIVATE ${columba_includes})
    target_compile_definitions(${name} PRIVATE "$<TARGET_PROPERTY:string_utils,INTERFACE_COMPILE_DEFINITIONS>")
    target_compile_options(${name} PRIVATE -O3 -fno-fast-math -fno-strict-aliasing)
    string_utils_silence_warnings(${name})
    target_link_libraries(${name} PRIVATE sdsl)
    set_target_properties(${name} PROPERTIES
      POSITION_INDEPENDENT_CODE ON LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bench/")
  endfunction()
  if(columba_bits_32)
    move_r_add_columba_plugin(columba-rlc-plugin-32
      "${columba_src}/bmove/bmove.cpp" "${columba_src}/bmove/moverepr.cpp")
    target_compile_definitions(columba-rlc-plugin-32 PRIVATE ALPHABET=5 RUN_LENGTH_COMPRESSION PHI_MOVE FMT_HEADER_ONLY THIRTY_TWO)
    move_r_add_columba_plugin(columba-fm-plugin-32
      "${columba_src}/fmindex/fmindex.cpp")
    target_compile_definitions(columba-fm-plugin-32 PRIVATE ALPHABET=5 FMT_HEADER_ONLY THIRTY_TWO)
  endif()
  if(columba_bits_64)
    move_r_add_columba_plugin(columba-rlc-plugin-64
      "${columba_src}/bmove/bmove.cpp" "${columba_src}/bmove/moverepr.cpp")
    target_compile_definitions(columba-rlc-plugin-64 PRIVATE ALPHABET=5 RUN_LENGTH_COMPRESSION PHI_MOVE FMT_HEADER_ONLY)
    move_r_add_columba_plugin(columba-fm-plugin-64
      "${columba_src}/fmindex/fmindex.cpp")
    target_compile_definitions(columba-fm-plugin-64 PRIVATE ALPHABET=5 FMT_HEADER_ONLY)
  endif()

  add_library(br_index INTERFACE)
  target_include_directories(br_index INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/external/br-index/src")
  target_link_libraries(br_index INTERFACE sdsl divsufsort divsufsort64)

  ############################# index-build CLI tools for the compared indexes #############################

  # br-index's own index builder
  add_executable(bri-build "${CMAKE_CURRENT_SOURCE_DIR}/external/br-index/src/bri-build.cpp")
  target_link_libraries(bri-build PRIVATE br_index malloc_count)
  set_target_properties(bri-build PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/")

  set(columba_build_sources
    "${columba_src}/buildindex.cpp" "${columba_src}/logger.cpp"
    "${columba_src}/parameters/buildparameters.cpp" "${columba_src}/seqfile.cpp")

  set(columba_bigbwt_dir "${CMAKE_CURRENT_SOURCE_DIR}/external/columba/external/Big-BWT")
  set(COLUMBA_BIGBWT_USABLE FALSE)
  if(EXISTS "${columba_bigbwt_dir}/utils.h")
    set(COLUMBA_BIGBWT_USABLE TRUE)
  endif()

  find_package(ZLIB REQUIRED)

  function(move_r_add_columba_build name rlc bits)
    add_executable(${name} ${columba_build_sources})
    target_include_directories(${name} SYSTEM PRIVATE ${columba_includes})
    target_compile_definitions(${name} PRIVATE "$<TARGET_PROPERTY:string_utils,INTERFACE_COMPILE_DEFINITIONS>")
    string_utils_silence_warnings(${name})
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
    add_executable(${name} "${CMAKE_CURRENT_SOURCE_DIR}/bench/move-rb-bench/indexes/columba-build-dispatch.cpp")
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
    list(APPEND columba_plugins columba-rlc-plugin-32 columba-fm-plugin-32)
  endif()
  if(columba_bits_64)
    list(APPEND columba_plugins columba-rlc-plugin-64 columba-fm-plugin-64)
  endif()

  ############################# sassy (online approximate matching) #############################
  # sassy (external/sassy) is a Rust library; cargo builds its C interface
  # (bench/move-rb-bench/indexes/sassy-capi) into bench/libsassy_capi.so, which move-rb-bench-apm loads.
  # It is only built if cargo and external/sassy exist.
  find_program(MOVE_R_CARGO cargo HINTS "$ENV{HOME}/.cargo/bin")
  set(sassy_capi_dir "${CMAKE_CURRENT_SOURCE_DIR}/bench/move-rb-bench/indexes/sassy-capi")
  if(MOVE_R_CARGO AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/sassy/Cargo.toml")
    set(sassy_capi_so "${CMAKE_BINARY_DIR}/bench/libsassy_capi.so")
    # sassy requires AVX2 or NEON at compile time (target-cpu=native)
    add_custom_command(OUTPUT "${sassy_capi_so}"
      COMMAND "${CMAKE_COMMAND}" -E env RUSTFLAGS=-Ctarget-cpu=native
              "${MOVE_R_CARGO}" build --release --manifest-path "${sassy_capi_dir}/Cargo.toml"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
              "${sassy_capi_dir}/target/release/libsassy_capi.so" "${sassy_capi_so}"
      DEPENDS "${sassy_capi_dir}/src/lib.rs" "${sassy_capi_dir}/Cargo.toml"
      COMMENT "Building sassy and its C ABI (cargo)" VERBATIM)
    add_custom_target(sassy-capi DEPENDS "${sassy_capi_so}")
  else()
    message(STATUS "sassy is not built (needs cargo and a checkout of external/sassy); "
                   "move-rb-bench-apm's --sassy will be unavailable")
  endif()

  move_r_add_tool(move-rb-bench-apm bench/move-rb-bench/move-rb-bench-apm.cpp
    "${CMAKE_BINARY_DIR}/bench/" br_index ${CMAKE_DL_LIBS})
  target_link_options(move-rb-bench-apm PRIVATE -rdynamic)
  add_dependencies(move-rb-bench-apm ${columba_plugins})
  if(TARGET sassy-capi)
    add_dependencies(move-rb-bench-apm sassy-capi)
  endif()
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
  # googletest and the test helpers (test-fuzz.hpp, test-inputs.hpp, test-references.hpp) are provided by string-utils
  include(GoogleTest)
  enable_testing()
  string_utils_add_components(tests)

  foreach(s move-r move-rb data-structures edit-distance-matrix)
    move_r_add_tool("test-${s}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test-${s}.cpp"
      "${CMAKE_BINARY_DIR}/tests/" string_utils_test)
    gtest_discover_tests("test-${s}")
  endforeach()

  foreach(s move-r move-rb data-structures edit-distance-matrix)
    set_target_properties("test-${s}" PROPERTIES INTERPROCEDURAL_OPTIMIZATION OFF)
  endforeach()
endif()

############################# silence the dependencies #############################

# disable the warnings of the external libraries
set(move_r_external_targets
  sparse_map rcomp r-index r-index-f OnlineRlbwt br_index
  grlbwt pfp-thresholds cdt bioparsers
  columba-rlc-plugin-32 columba-rlc-plugin-64 columba-fm-plugin-32 columba-fm-plugin-64
  columba-rlc-build-32 columba-rlc-build-64 columba-build-fm-32 columba-build-fm-64)

string_utils_mark_system(${move_r_external_targets})
string_utils_silence_warnings(${move_r_external_targets})
