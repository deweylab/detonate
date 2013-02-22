# LINUX
ifeq ($(shell uname), Linux)
  BOOST_INCLUDE = /ua/nathanae/downloads/boost/install/include
  BOOST_LIB     = /ua/nathanae/downloads/boost/install/lib
  BOOST_SUFFIX  = .a
  UNIT_TEST_DLL = libboost_unit_test_framework.so
  CC11          = /u/deweylab/sw/gcc-4.7.2/arch/x86_64-redhat-linux-gnu/bin/g++
endif

# MAC
ifeq ($(shell uname), Darwin)
  BOOST_INCLUDE = /opt/local/include
  BOOST_LIB     = /opt/local/lib
  BOOST_SUFFIX  = -mt.a
  UNIT_TEST_DLL = libboost_unit_test_framework-mt.dylib
endif

CC = /usr/bin/g++
#DEBUG = -g3 -fno-inline -O0
#CFLAGS = -O3 -W -Wall -Wextra $(DEBUG)
DEBUG =
CFLAGS = -O3 -fopenmp -W -Wall -Wextra $(DEBUG)
LFLAGS = -Wall $(DEBUG)
INCLUDE = -I$(BOOST_INCLUDE)
LIBS = $(BOOST_LIB)/libboost_program_options$(BOOST_SUFFIX) $(BOOST_LIB)/libboost_random$(BOOST_SUFFIX)
TEST_LIBS = $(BOOST_LIB)/$(UNIT_TEST_DLL) -Wl,-rpath,$(BOOST_LIB)/

all: summarize_kmer summarize

summarize_kmer: summarize_kmer.cpp summarize_kmer_meat.hh
	$(CC) $(CFLAGS) $(INCLUDE) summarize_kmer.cpp $(LIBS) -o summarize_kmer

summarize: summarize.cpp summarize_meat.hh
	#$(CC) $(CFLAGS) $(INCLUDE) summarize.cpp $(LIBS) -o summarize
	$(foreach gp, 1 2 3 4, \
	  $(foreach bp, 1 2 3 4, \
	    $(foreach np, 1 2, \
		  $(CC) $(CFLAGS) $(INCLUDE) -DGOOD_POLICY=$(gp) -DBETTER_POLICY=$(bp) -DN_POLICY=$(np) summarize.cpp $(LIBS) -o summarize_$(gp)_$(bp)_$(np);)))

test: test_lazycsv test_line_stream test_blast test_pslx test_psl test_pairset test_mask test_read_cluster_filter_alignments test_compute_alignment_stats
	./test_lazycsv
	./test_line_stream
	./test_blast
	./test_pslx
	./test_psl
	./test_pairset --show_progress
	./test_mask
	./test_read_cluster_filter_alignments
	./test_compute_alignment_stats

test_lazycsv: test_lazycsv.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_lazycsv.cpp $(LIBS) $(TEST_LIBS) -o test_lazycsv

test_line_stream: test_line_stream.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_line_stream.cpp $(LIBS) $(TEST_LIBS) -o test_line_stream

test_blast: test_blast.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_blast.cpp $(LIBS) $(TEST_LIBS) -o test_blast

test_pslx: test_pslx.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_pslx.cpp $(LIBS) $(TEST_LIBS) -o test_pslx

test_psl: test_psl.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_psl.cpp $(LIBS) $(TEST_LIBS) -o test_psl

test_pairset: test_pairset.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_pairset.cpp $(LIBS) $(TEST_LIBS) -o test_pairset

test_mask: test_mask.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_mask.cpp $(LIBS) $(TEST_LIBS) -o test_mask

test_read_cluster_filter_alignments: test_read_cluster_filter_alignments.cpp
	$(CC11) -std=c++11 $(CFLAGS) $(INCLUDE) test_read_cluster_filter_alignments.cpp $(LIBS) $(TEST_LIBS) -o test_read_cluster_filter_alignments

test_compute_alignment_stats: test_compute_alignment_stats.cpp
	$(CC11) -std=c++11 $(CFLAGS) $(INCLUDE) test_compute_alignment_stats.cpp $(LIBS) $(TEST_LIBS) -o test_compute_alignment_stats

# OLD

redux2: redux2.cpp
	$(CC) $(CFLAGS) $(INCLUDE) redux2.cpp $(LIBS) -o redux2

test_split: test_split.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_split.cpp $(LIBS) -o test_split

redux: redux.cpp
	$(CC) $(CFLAGS) $(INCLUDE) redux.cpp $(LIBS) -o redux

test_virtual: test_virtual.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_virtual.cpp $(LIBS) -o test_virtual

test_boost_icl: test_boost_icl.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_boost_icl.cpp $(LIBS) -o test_boost_icl

test_interval_for_blast: test_interval_for_blast.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_interval_for_blast.cpp $(LIBS) -o test_interval_for_blast

test_parse: test_parse.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_parse.cpp $(LIBS) -o test_parse

summarize_nucl_weighted_3: summarize_nucl_weighted_3.cpp psl.h
	$(CC) $(CFLAGS) $(INCLUDE) summarize_nucl_weighted_3.cpp $(LIBS) -o summarize_nucl_weighted_3

summarize_nucl_weighted_2: summarize_nucl_weighted_2.cpp psl.h
	$(CC) $(CFLAGS) $(INCLUDE) summarize_nucl_weighted_2.cpp $(LIBS) -o summarize_nucl_weighted_2
# DO NOT DELETE
