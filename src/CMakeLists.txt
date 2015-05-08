message(STATUS "Processing: main")

cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Enable verbose make output.") #TODO: remove
project(SGA C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
include(CheckLibraryExists)
include(CheckIncludeFileCXX)

# check that "SGA/sga.cpp" exists; if not, assume the source dir is wrong
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/SGA/sga.cpp")
    message(FATAL_ERROR "${PROJECT_SOURCE_DIR}: source files not found")
endif()

# set project-related variables that appear in config.h.in
#set(PACKAGE_VERSION "0.10.10")
set(PACKAGE_BUGREPORT "js18@sanger.ac.uk")
set(PACKAGE_URL "")

set(PACKAGE "${PROJECT_NAME}")
set(PACKAGE_NAME "${PROJECT_NAME}")
set(PACKAGE_TARNAME "${PROJECT_NAME}")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")

set(HEADER_SUBDIRS
    Algorithm Bigraph Concurrency GraphDiff Scaffold SGA SQG StringGraph SuffixTools Thirdparty Util
    tclap/include zstr/src
    CACHE INTERNAL "Subdirectories containing header files")
set(SUBDIRS
    BFTrim Dev SuffixTools Thirdparty Util
    CACHE INTERNAL "Subdirectories to descend into")

### Resolve external dependencies
#
# prefer static libraries
#set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so)

# bake-in RPATH to prevent library search problems
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
set(CMAKE_BUILD_WITH_INSTALL_RPATH true)

# check for OpenMP support
# => OPENMP_FOUND, OpenMP_C_FLAGS, OpenMP_CXX_FLAGS
if(NOT ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") AND NOT OpenMP_DETECTED)
    find_package(OpenMP)
    set(OpenMP_DETECTED TRUE CACHE BOOL "OpenMP detection was run")
    set(OPENMP_FOUND ${OPENMP_FOUND} CACHE BOOL "OpenMP found")
    set(OpenMP_C_FLAGS ${OpenMP_C_FLAGS} CACHE STRING "OpenMP C flags")
    set(OpenMP_CXX_FLAGS ${OpenMP_CXX_FLAGS} CACHE STRING "OpenMP CXX flags")
endif()

# find zlib
# => ZLIB_INCLUDE_DIRS, ZLIB_LIBRARIES
find_package(ZLIB REQUIRED)

# check for BamTools
# => BamTools_INCLUDE_DIR, BamTools_LIBRARY
find_package(BamTools REQUIRED)

# check for Google Sparsehash
# => SparseHash_INCLUDE_DIR
find_package(SparseHash REQUIRED)

# check for Boost 1.57
# => Boost_INCLUDE_DIRS, Boost_<COMPONENT>_LIBRARY
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND "${Boost_COMPILER}" STREQUAL "")
    message(WARNING
        "The cmake FindBoost package does not properly detect the compiler name when using CLang. Boost libraries built with 'versioned' layout (compiler name in library name, e.g. 'libboost_program_options-clang35-mt-1_57.so') might be incorrectly detected. To find such libraries, give the compiler name on the cmake command line, e.g. 'cmake -DBoost_COMPILER=\"-clang35\" [...]'. To check which libraries have been detected, use: 'grep libboost_program_options CMakeCache.txt'")
endif()
set(BOOST_ROOT "$ENV{BOOST_ROOT}" CACHE PATH "Path to Boost")
# Use static libraries to avoid LD_RUN_PATH issues
#set(Boost_USE_STATIC_LIBS true)
find_package(Boost 1.57.0 EXACT REQUIRED COMPONENTS iostreams program_options)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost 1.57.0 not found")
endif()

# check for Hoard
# => HOARD_FOUND, Hoard_LIBRARY
find_package(Hoard)

# check for JEMalloc
# => JEMALLOC_FOUND, JEMalloc_LIBRARY
find_package(JEMalloc)

# check for TCMalloc
# => TCMALLOC_FOUND, TCMalloc_LIBRARY
find_package(TCMalloc)

### Build configuration file
#
# check for clock_gettime()
# => HAVE_CLOCK_GETTIME
check_library_exists(rt clock_gettime "" HAVE_CLOCK_GETTIME)

# check for various hash map headers
# => HAVE_EXT_HASH_MAP, HAVE_GOOGLE_SPARSE_HASH_MAP, HAVE_TR1_UNORDERED_MAP, HAVE_UNORDERED_MAP
check_include_file_cxx(ext/hash_map HAVE_EXT_HASH_MAP)
check_include_file_cxx(google/sparse_hash_map HAVE_GOOGLE_SPARSE_HASH_MAP "-I${SparseHash_INCLUDE_DIR}")
check_include_file_cxx(tr1/unordered_map HAVE_TR1_UNORDERED_MAP)
check_include_file_cxx(unordered_map HAVE_UNORDERED_MAP)

# check for alloc libraries
# => HAVE_LIBHOARD, HAVE_LIBJEMALLOC, HAVE_LIBTCMALLOC
if(HOARD_FOUND)
    #check_library_exists(hoard malloc ${Hoard_LIBRARY} HAVE_LIBHOARD)
    set(HAVE_LIBHOARD 1)
    set(External_Malloc_LIBRARY ${Hoard_LIBRARY})
elseif(JEMALLOC_FOUND)
    #check_library_exists(jemalloc malloc ${JEMalloc_LIBRARY} HAVE_LIBJEMALLOC)
    set(HAVE_LIBJEMALLOC 1)
    set(External_Malloc_LIBRARY ${JEMalloc_LIBRARY})
elseif(TCMALLOC_FOUND)
    #check_library_exists(tcmalloc malloc ${TCMalloc_LIBRARY} HAVE_LIBTCMALLOC)
    set(HAVE_LIBTCMALLOC 1)
    set(External_Malloc_LIBRARY ${TCMalloc_LIBRARY})
endif()

# generate configuration file
configure_file(config.h.cmake.in config.h)

### Prepare default compilation flags
#
# headers
include_directories(SYSTEM
    ${ZLIB_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${BamTools_INCLUDE_DIR}
    ${SparseHash_INCLUDE_DIR})
include_directories(
    ${PROJECT_BINARY_DIR}
    ${HEADER_SUBDIRS})
get_directory_property(include_directories INCLUDE_DIRECTORIES)
message(STATUS "INCLUDE_DIRECTORIES='${include_directories}'")

### general compile flags
#
set(EXTRA_FLAGS "-Wall -Wextra -pedantic -pthread")

# compiler-specific flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -fmax-errors=1")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -ferror-limit=1")
endif()

# if OpenMP was found, add flags
#if(OPENMP_FOUND)
#    set(EXTRA_FLAGS "${EXTRA_FLAGS} ${OpenMP_C_FLAGS}")
#endif()

# consolidate flags
set(CMAKE_C_FLAGS "$ENV{CFLAGS} ${EXTRA_FLAGS}")
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} ${EXTRA_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS}'")

### build-specific flags
#
# compile flags for Debug build
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -ggdb -fno-inline -fno-eliminate-unused-debug-types")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
# compile flags for Test build
set(CMAKE_CXX_FLAGS_TEST "-O3 -g3")
set(CMAKE_C_FLAGS_TEST "${CMAKE_CXX_FLAGS_TEST}")
set(CMAKE_CXX_FLAGS_TEST_O2 "-O2 -g3")
set(CMAKE_C_FLAGS_TEST_O2 "${CMAKE_CXX_FLAGS_TEST_O2}")
set(CMAKE_CXX_FLAGS_TEST_O1 "-O1 -g3")
set(CMAKE_C_FLAGS_TEST_O1 "${CMAKE_CXX_FLAGS_TEST_O1}")
# compile flags for GProf build
set(CMAKE_CXX_FLAGS_GPROF "-O3 -g3 -pg")
set(CMAKE_C_FLAGS_GPROF "${CMAKE_CXX_FLAGS_GPROF}")
# compile flags for Release build
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DDISABLE_ASSERTS")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
# compile flags for GProfRel build
set(CMAKE_CXX_FLAGS_GPROFREL "-O3 -DNDEBUG -DDISABLE_ASSERTS -pg")
set(CMAKE_C_FLAGS_GPROFREL "${CMAKE_CXX_FLAGS_GPROF}")

# link flags
#set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=$ENV{LD_RUN_PATH} -Wl,--as-needed -static-libgcc -static-libstdc++")
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed")
message(STATUS "CMAKE_EXE_LINKER_FLAGS='${CMAKE_EXE_LINKER_FLAGS}'")

### Add target that updates version
#
add_custom_target(version ALL
    ${CMAKE_COMMAND} -DSRC_ROOT=${PROJECT_SOURCE_DIR} -DBIN_ROOT=${PROJECT_BINARY_DIR} -P ${PROJECT_SOURCE_DIR}/version.cmake)

### Enable testing targets
#
enable_testing()

### Descend into subdirectories
#
foreach(dir ${SUBDIRS})
    add_subdirectory(${dir})
endforeach()
