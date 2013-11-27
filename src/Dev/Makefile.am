# Note:
# This package defines default AM_C(XX)FLAGS and AM_LDFLAGS in configure.ac.
# Setting target-specific C(XX)FLAGS or LDFLAGS disables the use of the corresponding default.
AM_CXXFLAGS += -std=c++0x

AM_CPPFLAGS = \
	$(BOOST_CPPFLAGS) \
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

noinst_PROGRAMS = bloom-filter-tool test test_mac test_read_chunk_traversal
noinst_LTLIBRARIES = libmac.la libdev.la

bloom_filter_tool_LDFLAGS = $(AM_LDFLAGS) $(BOOST_IOSTREAMS_LDFLAGS)
bloom_filter_tool_LDADD = \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/Thirdparty/libthirdparty.a \
        libdev.la \
        $(BOOST_IOSTREAMS_LIBS)

bloom_filter_tool_SOURCES = bloom-filter-tool.cpp

test_CPPFLAGS = $(AM_CPPFLAGS) $(SPARSEHASH_CPPFLAGS)

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

test_SOURCES = test.cpp

libdev_la_SOURCES = kmer_gen.cpp

libmac_la_SOURCES = \
        Mutation.cpp \
        Read_Chunk.cpp \
        Read_Entry.cpp \
        Contig_Entry.cpp \
        Cigar.cpp \
        indent.cpp \
        Graph.cpp

test_mac_SOURCES = test_mac.cpp

test_mac_LDADD = libmac.la \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/Thirdparty/libthirdparty.a

test_read_chunk_traversal_SOURCES = test_read_chunk_traversal.cpp

test_read_chunk_traversal_LDADD = libmac.la \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/Thirdparty/libthirdparty.a
