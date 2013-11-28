bin_PROGRAMS = sga

sga_LDADD = \
        $(top_builddir)/Scaffold/libscaffold.a \
        $(top_builddir)/GraphDiff/libgraphdiff.a \
        $(top_builddir)/StringGraph/libstringgraph.a \
        $(top_builddir)/Concurrency/libconcurrency.a \
        $(top_builddir)/Algorithm/libalgorithm.a \
        $(top_builddir)/SuffixTools/libsuffixtools.a \
        $(top_builddir)/Bigraph/libbigraph.a \
        $(top_builddir)/Util/libutil.a \
        $(top_builddir)/SQG/libsqg.a \
        $(top_builddir)/Thirdparty/libthirdparty.a \
        $(BAMTOOLS_LIBS)

# Note:
# This package defines default AM_C(XX)FLAGS and AM_LDFLAGS in configure.ac.
# Setting target-specific C(XX)FLAGS or LDFLAGS disables the use of the corresponding default.
# Here, we override the default to enable pthread for sga only.
sga_LDFLAGS = $(AM_LDFLAGS) $(BAMTOOLS_LDFLAGS) -pthread

sga_SOURCES = sga.cpp \
        index.cpp index.h \
        overlap.cpp overlap.h \
        assemble.cpp assemble.h \
        correct.cpp correct.h \
        oview.cpp oview.h \
        preprocess.cpp preprocess.h \
        rmdup.cpp rmdup.h \
        merge.cpp merge.h \
        subgraph.cpp subgraph.h \
        scaffold.cpp scaffold.h \
        scaffold2fasta.cpp scaffold2fasta.h \
        connect.cpp connect.h \
        walk.cpp walk.h \
        filter.cpp filter.h \
        stats.cpp stats.h \
        fm-merge.cpp fm-merge.h \
        gmap.h gmap.cpp \
        filterBAM.h filterBAM.cpp \
        cluster.h cluster.cpp \
        gen-ssa.h gen-ssa.cpp \
        bwt2fa.h bwt2fa.cpp \
        graph-diff.h graph-diff.cpp \
        gapfill.h gapfill.cpp \
        preqc.h preqc.cpp \
        overlap-long.h overlap-long.cpp \
        variant-detectability.h variant-detectability.cpp \
        rewrite-evidence-bam.h rewrite-evidence-bam.cpp \
        haplotype-filter.h haplotype-filter.cpp \
        OverlapCommon.h OverlapCommon.cpp \
        SGACommon.h
