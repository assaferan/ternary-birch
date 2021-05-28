bin_PROGRAMS = birch birch_dbg
lib_LTLIBRARIES = libbirch.la

AM_CXXFLAGS = -O3 -DNDEBUG -Wall -Werror -std=c++11
DBG_AM_CXXFLAGS = -g -pg -Wall -Werror -std=c++11 -DDEBUG # -DDEBUG_LEVEL_FULL
AM_CXXFLAGS += -fvar-tracking-assignments-toggle
AM_CXXFLAGS += -fomit-frame-pointer
AM_CXXFLAGS += -funroll-all-loops
AM_LDFLAGS = -lm -lgmp -lgmpxx

SOURCES  = birch.h
SOURCES += birch_util.cpp
SOURCES += birch_util.h
SOURCES += Eigenvector.h
SOURCES += Fp.cpp
SOURCES += Fp.h
SOURCES += Genus.h
SOURCES += Genus.inl
SOURCES += HashMap.h
SOURCES += Isometry.h
SOURCES += Isometry.inl
SOURCES += IsometrySequence.h
SOURCES += Math.cpp
SOURCES += Math.h
SOURCES += Math.inl
SOURCES += Matrix.h
SOURCES += Matrix.inl
SOURCES += NeighborManager.h
SOURCES += NeighborManager.inl
SOURCES += ParseNipp.cpp
SOURCES += ParseNipp.h
SOURCES += Polynomial.h
SOURCES += Polynomial.inl
SOURCES += QuadForm.cpp
SOURCES += QuadForm.h
SOURCES += QuadForm.inl
SOURCES += Rational.h
SOURCES += Rational.inl
SOURCES += SetCover.cpp
SOURCES += SetCover.h
SOURCES += SquareMatrix.h
SOURCES += SquareMatrix.inl
SOURCES += Spinor.h

birch_dbg_SOURCES = birch.cpp $(SOURCES)
birch_dbg_CXXFLAGS = $(DBG_AM_CXXFLAGS)

birch_SOURCES = birch.cpp $(SOURCES)
birch_CXXFLAGS = $(AM_CXXFLAGS)

libbirch_la_SOURCES = $(SOURCES)
libbirch_la_CXXFLAGS = -shared -fPIC $(AM_CXXFLAGS)

birch_dbg_LDFLAGS = ${AM_LDFLAGS}

AM_LDFLAGS += -O3 -fPIC

birch_LDFLAGS = ${AM_LDFLAGS}
libbirch_la_LDFLAGS = ${AM_LDFLAGS}


