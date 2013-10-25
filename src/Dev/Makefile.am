noinst_PROGRAMS = bloom-filter-tool test test_mac
noinst_LIBRARIES = libmac.a libdev.a

bloom_filter_tool_CPPFLAGS = -std=c++0x -fopenmp \
        -I$(top_srcdir)/Util \
        -I$(top_srcdir)/Thirdparty

bloom_filter_tool_LDADD = \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/Thirdparty/libthirdparty.a \
        libdev.a \
        -lboost_iostreams

bloom_filter_tool_SOURCES = bloom-filter-tool.cpp

test_CPPFLAGS = -std=c++0x -fopenmp \
        -I$(top_srcdir)/Util \
        -I$(top_srcdir)/Bigraph \
        -I$(top_srcdir)/SuffixTools \
        -I$(top_srcdir)/StringGraph \
        -I$(top_srcdir)/Concurrency \
        -I$(top_srcdir)/Algorithm \
        -I$(top_srcdir)/SQG \
        -I$(top_srcdir)/Scaffold \
        -I$(top_srcdir)/GraphDiff \
        -I$(top_srcdir)/Thirdparty

test_LDADD = \
        $(top_builddir)/Scaffold/libscaffold.a \
        $(top_builddir)/GraphDiff/libgraphdiff.a \
        $(top_builddir)/StringGraph/libstringgraph.a \
        $(top_builddir)/Concurrency/libconcurrency.a \
        $(top_builddir)/Algorithm/libalgorithm.a \
        $(top_builddir)/SuffixTools/libsuffixtools.a \
        $(top_builddir)/Bigraph/libbigraph.a \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/SQG/libsqg.a \
        $(top_builddir)/Thirdparty/libthirdparty.a

test_SOURCES = \
        test.cpp

libdev_a_CPPFLAGS = -std=c++0x -fopenmp

libdev_a_SOURCES = kmer_gen.cpp

libmac_a_CPPFLAGS = -std=c++0x

libmac_a_SOURCES = \
        Mutation.cpp \
        Read_Chunk.cpp \
        Read_Entry.cpp \
        Contig_Entry.cpp \
        MAC.cpp

test_mac_CPPFLAGS = -std=c++0x

test_mac_SOURCES = test_mac.cpp

test_mac_LDADD = libmac.a
