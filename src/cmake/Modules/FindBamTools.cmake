# Find BamTools (https://github.com/pezmaster31/bamtools)
# Uses hint:
#   BAMTOOLS_ROOT
# Sets:
#   BAMTOOLS_FOUND
#   BamTools_INCLUDE_DIR
#   BamTools_LIBRARY
# Saves:
#   BAMTOOLS_ROOT
#   BamTools_INCLUDE_DIR_CACHED
#   BamTools_LIBRARY_CACHED


if(NOT BamTools_INCLUDE_DIR_CACHED OR NOT BamTools_LIBRARY_CACHED)
    set(BAMTOOLS_ROOT "$ENV{BAMTOOLS_ROOT}" CACHE PATH "Path to BamTools")

    find_path(BamTools_INCLUDE_DIR_CACHED api/api_global.h PATHS ${BAMTOOLS_ROOT}/include ${BAMTOOLS_ROOT}/include/bamtools)
    find_library(BamTools_LIBRARY_CACHED bamtools PATHS ${BAMTOOLS_ROOT}/lib ${BAMTOOLS_ROOT}/lib/bamtools)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(BamTools
        "BamTools library (https://github.com/pezmaster31/bamtools) not found. Specify location with -DBAMTOOLS_ROOT=<path>"
        BamTools_LIBRARY_CACHED BamTools_INCLUDE_DIR_CACHED)
    mark_as_advanced(BamTools_INCLUDE_DIR_CACHED BamTools_LIBRARY_CACHED)
else()
    message(STATUS "Using BamTools: ${BamTools_LIBRARY_CACHED}")
    set(BAMTOOLS_FOUND TRUE)
endif()

if(BAMTOOLS_FOUND)
    set(BamTools_INCLUDE_DIR ${BamTools_INCLUDE_DIR_CACHED})
    set(BamTools_LIBRARY ${BamTools_LIBRARY_CACHED})
endif()
