# Find SeqAn Library
# Uses hint:
#   SEQAN_ROOT
# Sets:
#   SEQAN_FOUND
#   SeqAn_INCLUDE_DIR
# Saves:
#   SEQAN_ROOT
#   SeqAn_INCLUDE_DIR_CACHED


if(NOT SeqAn_INCLUDE_DIR_CACHED)
    set(SEQAN_ROOT "$ENV{SEQAN_ROOT}" CACHE PATH "Path to SeqAn library.")

    find_path(SeqAn_INCLUDE_DIR_CACHED seqan/consensus.h ${SEQAN_ROOT}/include)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(SeqAn
        "SeqAn library (https://github.com/seqan/seqan) not found. Specify location with -DSEQAN_ROOT=<path>"
        SeqAn_INCLUDE_DIR_CACHED)
    mark_as_advanced(SeqAn_INCLUDE_DIR_CACHED)
else()
    message(STATUS "Using SeqAn: ${SeqAn_INCLUDE_DIR_CACHED}")
    set(SEQAN_FOUND TRUE)
endif()

if(SEQAN_FOUND)
    set(SeqAn_INCLUDE_DIR ${SeqAn_INCLUDE_DIR_CACHED})
endif()
