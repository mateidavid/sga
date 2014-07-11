message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(graphdiff STATIC
    DindelUtil.cpp
    DindelRealignWindow.cpp
    DindelHMM.cpp
    OverlapHaplotypeBuilder.cpp
    StringHaplotypeBuilder.cpp
    DeBruijnHaplotypeBuilder.cpp
    PairedDeBruijnHaplotypeBuilder.cpp
    HapgenUtil.cpp
    VCFTester.cpp
    GraphCompare.cpp
    OverlapExtractorWithCorrection.cpp
    DBGPathGuide.cpp
)
