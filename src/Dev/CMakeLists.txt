message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# adjust compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-ignored-qualifiers -Wno-unused-local-typedefs")
# order matters: local boost-intrusive must appear before general boost
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/interval-tree/include)
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/boost-intrusive/include) # <-first

add_library(mac STATIC
    globals.cpp
    global_assert.cpp
    Mutation_Chunk_Adapter.cpp
    Mutation.cpp
    Mutation_Cont.cpp
    Read_Chunk.cpp
    Read_Chunk_Cont.cpp
    Read_Entry.cpp
    Contig_Entry.cpp
    Cigar.cpp
    indent.cpp
    Graph.cpp
)
#set_target_properties(mac PROPERTIES LINK_FLAGS "${BOOST_LDFLAGS}")

# workaround library&executable having the same name
add_executable(mac-bin mac.cpp)
set_target_properties(mac-bin PROPERTIES OUTPUT_NAME mac)
add_dependencies(mac-bin version)
#set_target_properties(mac-bin PROPERTIES LINK_FLAGS "${BOOST_LDFLAGS} -lboost_iostreams -lboost_program_options")
#target_link_libraries(mac-bin mac ${CMAKE_BINARY_DIR}/Util/libutil.a ${CMAKE_BINARY_DIR}/Thirdparty/libthirdparty.a)
target_link_libraries(mac-bin
    mac
    util
    thirdparty
    ${Boost_IOSTREAMS_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${ZLIB_LIBRARIES}
)

add_executable(test-itree-and-factory
    test-itree-and-factory.cpp
    globals.cpp
    global_assert.cpp
)
target_link_libraries(test-itree-and-factory
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

file(COPY factory.gdb mac-gdb.gdb DESTINATION .)
