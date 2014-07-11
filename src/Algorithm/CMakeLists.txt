message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(algorithm STATIC
    OverlapAlgorithm.cpp
    SearchSeed.cpp
    OverlapBlock.cpp
    SearchHistory.cpp
    ErrorCorrectProcess.cpp
    QCProcess.cpp
    OverlapTools.cpp
    DPAlignment.cpp
    ConnectProcess.cpp
    StringGraphGenerator.cpp
    FMMergeProcess.cpp
    StatsProcess.cpp
    ClusterProcess.cpp
    ReadCluster.cpp
    LRAlignment.cpp
    StringThreader.cpp
    ExtensionDP.cpp
    HaplotypeBuilder.cpp
    GapFillProcess.cpp
    VariationBuilderCommon.cpp
    KmerOverlaps.cpp
)
