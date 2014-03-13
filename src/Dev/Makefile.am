# Note:
# This package defines default AM_C(XX)FLAGS and AM_LDFLAGS in configure.ac.
# Setting target-specific C(XX)FLAGS or LDFLAGS disables the use of the corresponding default.
AM_CXXFLAGS += -std=c++0x

noinst_PROGRAMS = bloom-filter-tool mac test-itree-and-factory
noinst_LTLIBRARIES = libmac.la libdev.la

bloom_filter_tool_LDFLAGS = $(AM_LDFLAGS) $(BOOST_LDFLAGS)
bloom_filter_tool_LDADD = \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/Thirdparty/libthirdparty.a \
        libdev.la \
        $(BOOST_IOSTREAMS_LIBS)

bloom_filter_tool_SOURCES = bloom-filter-tool.cpp

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

#test_SOURCES = test.cpp

libdev_la_SOURCES = kmer_gen.cpp

libmac_la_SOURCES = \
	globals.cpp \
	common.cpp \
        Mutation.cpp \
        Read_Chunk.cpp \
        Read_Entry.cpp \
        Contig_Entry.cpp \
        Cigar.cpp \
        indent.cpp \
        Graph.cpp

mac_SOURCES = mac.cpp
mac_LDFLAGS = $(AM_LDFLAGS) $(BOOST_LDFLAGS)
mac_LDADD = libmac.la \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/Thirdparty/libthirdparty.a \
	$(BOOST_IOSTREAMS_LIBS)

#test_read_chunk_traversal_SOURCES = test_read_chunk_traversal.cpp
#test_read_chunk_traversal_LDADD = libmac.la \
#	$(top_builddir)/Util/libutil.a \
#	$(top_builddir)/Thirdparty/libthirdparty.a

test_itree_and_factory_SOURCES = test-itree-and-factory.cpp common.cpp globals.cpp
test_itree_and_factory_CPPFLAGS = -I$(top_srcdir)/boost-intrusive/include -I$(top_srcdir)/interval-tree/include $(AM_CPPFLAGS)
test_itree_and_factory_CXXFLAGS = $(AM_CXXFLAGS) -Wno-ignored-qualifiers
test_itree_and_factory_LDFLAGS = $(AM_LDFLAGS) $(BOOST_LDFLAGS)
test_itree_and_factory_LDADD = $(BOOST_PROGRAM_OPTIONS_LIBS)
