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
  CC11          = g++
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

.PHONY: all
all: summarize summarize_kmer

summarize_jobs := $(foreach gp, 1 2 3 4, $(foreach bp, 1 2 3 4, $(foreach np, 1 2, summarize_${gp}_${bp}_${np})))
gp = $(word 1,$(subst _, ,$*))
bp = $(word 2,$(subst _, ,$*))
np = $(word 3,$(subst _, ,$*))

.PHONY: summarize
summarize: ${summarize_jobs}
${summarize_jobs}: summarize_%: summarize.cpp summarize_meat.hh
	$(CC) $(CFLAGS) $(INCLUDE) -DGOOD_POLICY=$(gp) -DBETTER_POLICY=$(bp) -DN_POLICY=$(np) summarize.cpp $(LIBS) -o summarize_$(gp)_$(bp)_$(np)

summarize_kmer: summarize_kmer.cpp summarize_kmer_meat.hh
	$(CC) $(CFLAGS) $(INCLUDE) summarize_kmer.cpp $(LIBS) -o summarize_kmer

all_tests := test_lazycsv test_line_stream test_blast test_pslx test_psl test_pairset test_mask test_read_cluster_filter_alignments test_compute_alignment_stats test_alignment_segment

test: ${all_tests}
	./test_lazycsv
	./test_line_stream
	./test_blast
	./test_pslx
	./test_psl
	./test_pairset --show_progress
	./test_mask
	./test_read_cluster_filter_alignments
	./test_compute_alignment_stats
	./test_alignment_segment

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

test_alignment_segment: test_alignment_segment.cpp alignment_segment.hh
	$(CC11) -std=c++11 $(CFLAGS) $(INCLUDE) test_alignment_segment.cpp $(LIBS) $(TEST_LIBS) -o test_alignment_segment

.PHONY:
clean:
	-rm -f ${summarize_jobs} summarize_kmer ${all_tests}
