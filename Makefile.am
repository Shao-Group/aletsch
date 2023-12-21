AUTOMAKE_OPTIONS = foreign
EXTRA_DIST = LICENSE
SUBDIRS = util graph gtf rnacore scallop bridge meta

bin_PROGRAMS = aletsch

UTIL_INCLUDE = $(top_srcdir)/util
GTF_INCLUDE = $(top_srcdir)/gtf
GRAPH_INCLUDE = $(top_srcdir)/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore
SCALLOP_INCLUDE = $(top_srcdir)/scallop
BRIDGE_INCLUDE = $(top_srcdir)/bridge
META_INCLUDE = $(top_srcdir)/meta

UTIL_LIB = $(top_builddir)/util
GTF_LIB = $(top_builddir)/gtf
GRAPH_LIB = $(top_builddir)/graph
RNACORE_LIB = $(top_builddir)/rnacore
SCALLOP_LIB = $(top_builddir)/scallop
BRIDGE_LIB = $(top_builddir)/bridge
META_LIB = $(top_builddir)/meta

aletsch_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(SCALLOP_INCLUDE) -I$(RNACORE_INCLUDE) -I$(BRIDGE_INCLUDE) -I$(META_INCLUDE)
aletsch_LDFLAGS = -pthread -L$(GTF_LIB) -L$(GRAPH_LIB) -L$(UTIL_LIB) -L$(SCALLOP_LIB) -L$(RNACORE_LIB) -L$(BRIDGE_LIB) -L$(META_LIB)
aletsch_LDADD = -lmeta -lscallop -lbridge -lrnacore -lgtf -lgraph -lutil 
aletsch_SOURCES = aletsch.cc
