# ------------------
# - Configuration. -
# ------------------

# You should edit the following variables as needed for your system.
#
# Notes:
#
# - CC11 is only needed to build the tests. You can ignore it if you just want
#   to build the main REF-EVAL package.
#
# - OMP should be -fopenmp if OpenMP is suppported by your compiler, and you
#   want to be able to run in multithreaded mode. OpenMP is supported by GCC
#   and by the latest version of Clang, but not (as of this writing) by the
#   version of Clang distributed by Apple.

ifeq ($(shell uname), Linux)
  CC   = g++
  CC11 = $(CC) -std=c++11 
  OMP  = -fopenmp
endif

ifeq ($(shell uname), Darwin)
  CC   = clang++
  CC11 = $(CC) -std=c++11 
  OMP  =
endif

# ------------------------
# - Build specification. -
# ------------------------

# You shouldn't normally need to edit anything below here.

#DEBUG = -g3 -fno-inline -O0 -Wall -Wextra 
DEBUG =
CFLAGS = -W -g -O3 $(DEBUG)
BOOST_INCLUDE = -isystem boost
BOOST_LIB     = boost/stage/lib/libboost_program_options.a boost/stage/lib/libboost_random.a
LEMON_INCLUDE = -Ilemon/build -Ilemon/lemon-main-473c71baff72
LEMON_LIB     = lemon/build/lemon/libemon.a -lpthread
CITY_INCLUDE  = -Icity/install/include
CITY_LIB      = city/install/lib/libcityhash.a
SH_INCLUDE 	  = -Isparsehash/include
DL_INCLUDE 	  = -Ideweylab
INCLUDE = $(BOOST_INCLUDE) $(LEMON_INCLUDE) $(CITY_INCLUDE) $(SH_INCLUDE) $(DL_INCLUDE)
LIBS    = $(BOOST_LIB) $(LEMON_LIB) $(CITY_LIB)
TEST_LIBS = -lboost_unit_test_framework

.PHONY: all
all: ref-eval

boost/finished:
	@echo 
	@echo -----------------------------------------------
	@echo - Building a subset of the Boost C++ library. -
	@echo -----------------------------------------------
	@echo 
	cd boost && $(MAKE)

lemon/finished:
	@echo 
	@echo -------------------------------------
	@echo - Building the Lemon graph library. -
	@echo -------------------------------------
	@echo 
	cd lemon && $(MAKE)

city/finished:
	@echo 
	@echo ----------------------------------
	@echo - Building the CityHash library. -
	@echo ----------------------------------
	@echo 
	cd city && $(MAKE)

sparsehash/finished:
	@echo 
	@echo ------------------------------------
	@echo - Building the SparseHash library. -
	@echo ------------------------------------
	@echo 
	cd sparsehash && $(MAKE)

ref-eval: ref-eval.cpp boost/finished lemon/finished city/finished sparsehash/finished
	@echo 
	@echo -----------------------------
	@echo - Building REF-EVAL itself. -
	@echo -----------------------------
	@echo 
	$(CC) $(OMP) $(CFLAGS) $(INCLUDE) ref-eval.cpp $(LIBS) -o ref-eval

all_tests := test_lazycsv test_line_stream test_blast test_psl test_pairset test_mask test_read_cluster_filter_alignments test_compute_alignment_stats test_alignment_segment test_re_matched

.PHONY: test
test: ${all_tests} boost/finished lemon/finished
	./test_lazycsv
	./test_line_stream
	./test_blast
	./test_psl
	./test_pairset --show_progress
	./test_mask
	./test_read_cluster_filter_alignments
	./test_compute_alignment_stats
	./test_alignment_segment
	./test_re_matched

test_lazycsv: test_lazycsv.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_lazycsv.cpp $(LIBS) $(TEST_LIBS) -o test_lazycsv

test_line_stream: test_line_stream.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_line_stream.cpp $(LIBS) $(TEST_LIBS) -o test_line_stream

test_blast: test_blast.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_blast.cpp $(LIBS) $(TEST_LIBS) -o test_blast

test_psl: test_psl.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_psl.cpp $(LIBS) $(TEST_LIBS) -o test_psl

test_pairset: test_pairset.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_pairset.cpp $(LIBS) $(TEST_LIBS) -o test_pairset

test_mask: test_mask.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_mask.cpp $(LIBS) $(TEST_LIBS) -o test_mask

test_read_cluster_filter_alignments: test_read_cluster_filter_alignments.cpp
	$(CC11) $(CFLAGS) $(INCLUDE) test_read_cluster_filter_alignments.cpp $(LIBS) $(TEST_LIBS) -o test_read_cluster_filter_alignments

test_compute_alignment_stats: test_compute_alignment_stats.cpp
	$(CC11) $(CFLAGS) $(INCLUDE) test_compute_alignment_stats.cpp $(LIBS) $(TEST_LIBS) -o test_compute_alignment_stats

test_alignment_segment: test_alignment_segment.cpp alignment_segment.hh
	$(CC11) $(CFLAGS) $(INCLUDE) test_alignment_segment.cpp $(LIBS) $(TEST_LIBS) -o test_alignment_segment

test_re_matched: test_re_matched.cpp re_matched_meat.hh
	$(CC11) $(CFLAGS) $(INCLUDE) test_re_matched.cpp $(LIBS) $(TEST_LIBS) -o test_re_matched

.PHONY:
clean:
	-rm -f ref-eval ${all_tests}
	-cd boost && $(MAKE) clean
	-cd lemon && $(MAKE) clean
	-cd city && $(MAKE) clean
	-cd sparsehash && $(MAKE) clean
