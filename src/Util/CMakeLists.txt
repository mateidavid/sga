message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(util STATIC
    Util.cpp
    Alphabet.cpp
    Contig.cpp
    ReadTable.cpp
    ReadInfoTable.cpp
    SeqReader.cpp
    DNAString.cpp
    Match.cpp
    Pileup.cpp
    BitChar.cpp
    Interval.cpp
    SeqCoord.cpp
    MultiOverlap.cpp
    QualityVector.cpp
    SGAStats.cpp
    Quality.cpp
    PrimerScreen.cpp
    BitVector.cpp
    CorrectionThresholds.cpp
    KmerDistribution.cpp
    ClusterReader.cpp
    MultiAlignment.cpp
    StdAlnTools.cpp
    VCFUtil.cpp
    QualityTable.cpp
    BloomFilter.cpp
)
