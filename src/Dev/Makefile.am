AM_CXXFLAGS = $(BOOST_CPPFLAGS)

noinst_PROGRAMS = bloom-filter-tool test test_mac test_read_chunk_traversal
noinst_LTLIBRARIES = libmac.la libdev.la

bloom_filter_tool_CPPFLAGS = -std=c++0x -fopenmp \
        -I$(top_srcdir)/Util \
        -I$(top_srcdir)/Thirdparty

bloom_filter_tool_LDADD = \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/Thirdparty/libthirdparty.a \
        libdev.la \
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

libdev_la_CPPFLAGS = -std=c++0x -fopenmp

libdev_la_SOURCES = kmer_gen.cpp

libmac_la_CPPFLAGS = -std=c++0x \
        -I$(top_srcdir)/Bigraph \
        -I$(top_srcdir)/Thirdparty \
        -I$(top_srcdir)/Util

libmac_la_SOURCES = \
        Mutation.cpp \
        Read_Chunk.cpp \
        Read_Entry.cpp \
        Contig_Entry.cpp \
        Cigar.cpp \
        indent.cpp \
        Graph.cpp

test_mac_CPPFLAGS = -std=c++0x

test_mac_SOURCES = test_mac.cpp

test_mac_LDADD = libmac.la \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/Thirdparty/libthirdparty.a

test_read_chunk_traversal_CPPFLAGS = -std=c++0x

test_read_chunk_traversal_SOURCES = test_read_chunk_traversal.cpp

test_read_chunk_traversal_LDADD = libmac.la \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/Thirdparty/libthirdparty.a
