ifeq ($(shell uname), Linux)
  CC11 = /u/deweylab/sw/gcc-4.7.2/arch/x86_64-redhat-linux-gnu/bin/g++
  CC = g++ -W
  OMP = -fopenmp
endif
ifeq ($(shell uname), Darwin)
  CC11 = clang++ -W -Wno-unused-parameter 
  CC   = clang++ -W -Wno-unused-parameter
  OMP  =
endif

#DEBUG = -g3 -fno-inline -O0 -Wall -Wextra 
DEBUG =
CFLAGS = -g -O3 $(DEBUG)
SH_INCLUDE 	  = -Isparsehash
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
all: summarize summarize_axolotl summarize_matched summarize_oomatched summarize_aligned_kmer summarize_kmer summarize_wkr summarize_kmerpair summarize_multikmer

force_look:
	true

boost/finished:
	cd boost && ./bootstrap.sh && ./b2 && touch finished

lemon/finished:
	cd lemon && make

city/finished:
	cd city && make

summarize_jobs := $(foreach gp, 1 2 3 4 5 6, $(foreach bp, 1 2 3 4 5, $(foreach np, 1 2, $(foreach mpi, 80 95, summarize_${gp}_${bp}_${np}_${mpi}))))
gp = $(word 1,$(subst _, ,$*))
bp = $(word 2,$(subst _, ,$*))
np = $(word 3,$(subst _, ,$*))
mpi = $(word 4,$(subst _, ,$*))

ref-eval: ref-eval.cpp boost/finished lemon/finished city/finished
	$(CC) $(OMP) $(CFLAGS) $(INCLUDE) ref-eval.cpp $(LIBS) -o ref-eval

.PHONY: summarize
summarize: ${summarize_jobs} boost/finished lemon/finished
${summarize_jobs}: summarize_%: summarize.cpp summarize_meat.hh
	$(CC) $(OMP) $(CFLAGS) $(INCLUDE) -DGOOD_POLICY=$(gp) -DBETTER_POLICY=$(bp) -DN_POLICY=$(np) -DMIN_PCT_ID=$(mpi) summarize.cpp $(LIBS) -o summarize_$(gp)_$(bp)_$(np)_$(mpi)

summarize_axolotl: summarize_axolotl.cpp summarize_axolotl_meat.hh boost/finished lemon/finished
	$(CC) $(OMP) $(CFLAGS) $(INCLUDE) summarize_axolotl.cpp $(LIBS) -o summarize_axolotl

summarize_matched: summarize_matched.cpp summarize_matched_meat.hh boost/finished lemon/finished
	$(CC) $(OMP) $(CFLAGS) $(INCLUDE) summarize_matched.cpp $(LIBS) -o summarize_matched

summarize_oomatched: summarize_oomatched.cpp summarize_oomatched_meat.hh boost/finished lemon/finished
	$(CC) $(CFLAGS) $(INCLUDE) summarize_oomatched.cpp $(LIBS) -o summarize_oomatched

summarize_kmer: summarize_kmer.cpp summarize_kmer_meat.hh boost/finished lemon/finished
	$(CC) $(CFLAGS) $(INCLUDE) summarize_kmer.cpp $(LIBS) -o summarize_kmer

summarize_wkr: summarize_wkr.cpp summarize_wkr_meat.hh boost/finished lemon/finished
	$(CC) $(CFLAGS) $(INCLUDE) summarize_wkr.cpp $(LIBS) -o summarize_wkr

summarize_aligned_kmer: summarize_aligned_kmer.cpp summarize_aligned_kmer_meat.hh boost/finished lemon/finished
	$(CC) $(CFLAGS) $(INCLUDE) summarize_aligned_kmer.cpp $(LIBS) -o summarize_aligned_kmer

summarize_multikmer: summarize_multikmer.cpp summarize_multikmer_meat.hh boost/finished lemon/finished
	$(CC) $(OMP) $(CFLAGS) $(INCLUDE) summarize_multikmer.cpp $(LIBS) -o summarize_multikmer

summarize_kmerpair: summarize_kmerpair.cpp summarize_kmerpair_meat.hh boost/finished lemon/finished
	$(CC) $(CFLAGS) $(INCLUDE) summarize_kmerpair.cpp $(LIBS) -o summarize_kmerpair

all_tests := test_lazycsv test_line_stream test_blast test_psl test_pairset test_mask test_read_cluster_filter_alignments test_compute_alignment_stats test_alignment_segment test_summarize_matched

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
	./test_summarize_matched

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
	$(CC11) -std=c++11 $(CFLAGS) $(INCLUDE) test_read_cluster_filter_alignments.cpp $(LIBS) $(TEST_LIBS) -o test_read_cluster_filter_alignments

test_compute_alignment_stats: test_compute_alignment_stats.cpp
	$(CC11) -std=c++11 $(CFLAGS) $(INCLUDE) test_compute_alignment_stats.cpp $(LIBS) $(TEST_LIBS) -o test_compute_alignment_stats

test_alignment_segment: test_alignment_segment.cpp alignment_segment.hh
	$(CC11) -std=c++11 $(CFLAGS) $(INCLUDE) test_alignment_segment.cpp $(LIBS) $(TEST_LIBS) -o test_alignment_segment

test_summarize_matched: test_summarize_matched.cpp summarize_matched_meat.hh
	$(CC11) -std=c++11 $(CFLAGS) $(INCLUDE) test_summarize_matched.cpp $(LIBS) $(TEST_LIBS) -o test_summarize_matched

.PHONY:
clean:
	-rm -f ${summarize_jobs} summarize_axolotl summarize_matched summarize_oomatched summarize_aligned_kmer summarize_kmer summarize_wkr summarize_kmerpair summarize_multikmer ${all_tests}
	cd lemon && make clean
	cd city && make clean
