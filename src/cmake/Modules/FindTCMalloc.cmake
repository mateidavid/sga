# Find TCMalloc
# Uses hint:
#   TCMALLOC_ROOT
# Sets:
#   TCMALLOC_FOUND
#   TCMalloc_LIBRARY
# Saves:
#   TCMALLOC_ROOT
#   TCMalloc_LIBRARY_CACHED


if(NOT TCMalloc_LIBRARY_CACHED)
    set(TCMALLOC_ROOT "$ENV{TCMALLOC_ROOT}" CACHE PATH "Path to TCMalloc")

    find_library(TCMalloc_LIBRARY_CACHED NAMES tcmalloc tcmalloc_minimal PATHS ${TCMALLOC_ROOT} NO_DEFAULT_PATH)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(TCMalloc
        "TCMalloc library (http://goog-perftools.sourceforge.net/doc/tcmalloc.html) not found. Specify location with -DTCMALLOC_ROOT=<path>"
        TCMalloc_LIBRARY_CACHED)
    mark_as_advanced(TCMalloc_LIBRARY_CACHED)
else()
    message(STATUS "Using TCMalloc: ${TCMalloc_LIBRARY_CACHED}")
    set(TCMALLOC_FOUND TRUE)
endif()

if(TCMALLOC_FOUND)
    set(TCMalloc_LIBRARY ${TCMalloc_LIBRARY_CACHED})
endif()
