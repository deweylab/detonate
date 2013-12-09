ifeq ($(shell uname), Linux)
  CC = g++ -W
  CC11 = g++ -std=c++11 
  OMP = -fopenmp
  HASH_FUN_H = "<tr1/functional>"
  HASH_NAMESPACE = "std::tr1"
endif
ifeq ($(shell uname), Darwin)
  CC11 = clang++ -W -Wno-unused-parameter 
  CC   = clang++ -W -Wno-unused-parameter
  OMP  =
  HASH_FUN_H = "<functional>"
  HASH_NAMESPACE = "std"
endif

#DEBUG = -g3 -fno-inline -O0 -Wall -Wextra 
DEBUG =
CFLAGS = -g -O3 $(DEBUG)
SH_INCLUDE 	  = -Isparsehash -DHASH_FUN_H=$(HASH_FUN_H) -DHASH_NAMESPACE=$(HASH_NAMESPACE)
DL_INCLUDE 	  = -Ideweylab
BOOST_INCLUDE = -Iboost
BOOST_LIB     = boost/stage/lib/libboost_program_options.a boost/stage/lib/libboost_random.a
LEMON_INCLUDE = -Ilemon/build -Ilemon/lemon-main-473c71baff72
LEMON_LIB     = lemon/build/lemon/libemon.a -lpthread
CITY_INCLUDE  = -Icity/install/include
CITY_LIB      = city/install/lib/libcityhash.a
INCLUDE = $(SH_INCLUDE) $(DL_INCLUDE) $(BOOST_INCLUDE) $(LEMON_INCLUDE) $(CITY_INCLUDE)
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

ref-eval: ref-eval.cpp boost/finished lemon/finished city/finished
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
