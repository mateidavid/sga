# Find JEMalloc
# Uses hint:
#   JEMALLOC_ROOT
# Sets:
#   JEMALLOC_FOUND
#   JEMalloc_LIBRARY
# Saves:
#   JEMALLOC_ROOT
#   JEMalloc_LIBRARY_CACHED


if(NOT JEMalloc_LIBRARY_CACHED)
    set(JEMALLOC_ROOT "$ENV{JEMALLOC_ROOT}" CACHE PATH "Path to JEMalloc")

    find_library(JEMalloc_LIBRARY_CACHED jemalloc PATHS ${JEMALLOC_ROOT} NO_DEFAULT_PATH)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(JEMalloc
        "JEMalloc library (http://www.canonware.com/jemalloc) not found. Specify location with -DJEMALLOC_ROOT=<path>"
        JEMalloc_LIBRARY_CACHED)
    mark_as_advanced(JEMalloc_LIBRARY_CACHED)
else()
    message(STATUS "Using JEMalloc: ${JEMalloc_LIBRARY_CACHED}")
    set(JEMALLOC_FOUND TRUE)
endif()

if(JEMALLOC_FOUND)
    set(JEMalloc_LIBRARY ${JEMalloc_LIBRARY_CACHED})
endif()
