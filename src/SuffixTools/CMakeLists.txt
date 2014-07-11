message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(suffixtools STATIC
    STCommon.cpp 
    Occurrence.cpp 
    SuffixArray.cpp 
    SuffixCompare.cpp 
    InverseSuffixArray.cpp 
    SACAInducedCopying.cpp
    BWTAlgorithms.cpp
    BWTDiskConstruction.cpp
    BWTReader.cpp
    BWTWriter.cpp
    SAReader.cpp
    SAWriter.cpp
    GapArray.cpp
    RankProcess.cpp
    SBWT.cpp
    RLBWT.cpp
    BWTReader.cpp
    BWTWriter.cpp
    BWTWriterBinary.cpp
    BWTReaderBinary.cpp
    BWTWriterAscii.cpp
    BWTReaderAscii.cpp
    BWTIntervalCache.cpp
    QuickBWT.cpp
    SampledSuffixArray.cpp
    BWTCABauerCoxRosone.cpp
    BWTCARopebwt.cpp
    PopulationIndex.cpp
)
