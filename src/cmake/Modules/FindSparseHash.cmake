# Find Google SparseHash
# Uses hint:
#   SPARSEHASH_ROOT
# Sets:
#   SPARSEHASH_FOUND
#   SparseHash_INCLUDE_DIR
# Saves:
#   SPARSEHASH_ROOT
#   SparseHash_INCLUDE_DIR_CACHED


if(NOT SparseHash_INCLUDE_DIR_CACHED)
    set(SPARSEHASH_ROOT "$ENV{SPARSEHASH_ROOT}" CACHE PATH "Path to Google SparseHash library.")

    find_path(SparseHash_INCLUDE_DIR_CACHED google/sparse_hash_map ${SPARSEHASH_ROOT}/include)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(SparseHash
        "Google SparseHash library (http://code.google.com/p/google-sparsehash) not found. Specify location with -DSPARSEHASH_ROOT=<path>"
        SparseHash_INCLUDE_DIR_CACHED)
    mark_as_advanced(SparseHash_INCLUDE_DIR_CACHED)
else()
    message(STATUS "Using SparseHash: ${SparseHash_INCLUDE_DIR_CACHED}")
    set(SPARSEHASH_FOUND TRUE)
endif()

if(SPARSEHASH_FOUND)
    set(SparseHash_INCLUDE_DIR ${SparseHash_INCLUDE_DIR_CACHED})
endif()
