noinst_LIBRARIES = libthirdparty.a

libthirdparty_a_SOURCES = \
        gzstream.h gzstream.C \
        stdaln.h stdaln.c \
        overlapper.h overlapper.cpp \
        multiple_alignment.h multiple_alignment.cpp \
        bcr.h bcr.c \
        MurmurHash3.h MurmurHash3.cpp \
        count_min_sketch.h count_min_sketch.cpp \
        rapidjson/allocators.h \
        rapidjson/document.h \
        rapidjson/encodedstream.h \
        rapidjson/encodings.h \
        rapidjson/filereadstream.h \
        rapidjson/filestream.h \
        rapidjson/filewritestream.h \
        rapidjson/prettywriter.h \
        rapidjson/rapidjson.h \
        rapidjson/reader.h \
        rapidjson/stringbuffer.h \
        rapidjson/writer.h \
        rapidjson/internal/pow10.h \
        rapidjson/internal/stack.h \
        rapidjson/internal/strfunc.h

# Note:
# This package defines default AM_C(XX)FLAGS and AM_LDFLAGS in configure.ac.
# Setting target-specific C(XX)FLAGS or LDFLAGS disables the use of the corresponding default.
# Here, we override the default to disable -Werror for third-party libraries only.
libthirdparty_a_CFLAGS = $(AM_CFLAGS) -Wno-error
libthirdparty_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-error
