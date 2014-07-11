# Find Hoard
# Uses hint:
#   HOARD_ROOT
# Sets:
#   HOARD_FOUND
#   Hoard_LIBRARY
# Saves:
#   HOARD_ROOT
#   Hoard_LIBRARY_CACHED


if(NOT Hoard_LIBRARY_CACHED)
    set(HOARD_ROOT "$ENV{HOARD_ROOT}" CACHE PATH "Path to Hoard")

    find_library(Hoard_LIBRARY_CACHED hoard PATHS ${HOARD_ROOT} NO_DEFAULT_PATH)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Hoard
        "Hoard library (http://www.hoard.org/) not found. Specify location with -DHOARD_ROOT=<path>"
        Hoard_LIBRARY_CACHED)
    mark_as_advanced(Hoard_LIBRARY_CACHED)
else()
    message(STATUS "Using Hoard: ${Hoard_LIBRARY_CACHED}")
    set(HOARD_FOUND TRUE)
endif()

if(HOARD_FOUND)
    set(Hoard_LIBRARY ${Hoard_LIBRARY_CACHED})
endif()
