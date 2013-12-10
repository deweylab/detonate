# ------------------
# - Configuration. -
# ------------------

# You should edit the following variables as needed for your system. It is very
# likely that you don't need to edit anything if you just want to build the
# main ref-eval executable on Linux or Mac. In that case, just type "make" and
# wait for everything to compile.
#
# Notes:
#
# - CXX11 is only needed to build the tests. You can ignore it if you just want
#   to build the main REF-EVAL package. On some systems, the default GCC is
#   somewhat old (e.g., GCC 4.4.6 on RHEL 6.3), and might not support C++11. In
#   that case, either just run the tests somewhere else or install a newer
#   version of GCC and update CXX11 accordingly.
#
# - OMP should be -fopenmp if OpenMP is suppported by your compiler, and you
#   want to be able to run in multithreaded mode. OpenMP is supported by GCC,
#   but not (as of this writing) by the version of Clang distributed by Apple.
#   If you want to use -fopenmp on a Mac, it is easy to install GCC from
#   MacPorts (http://www.macports.org/). If you do so, you'll want to update
#   CXX below to use MacPorts' g++ instead of clang++.

ifeq ($(shell uname), Linux)
  CXX   = g++
  CXX11 = $(CXX) -std=c++11 
  OMP   = -fopenmp
endif

ifeq ($(shell uname), Darwin)
  CXX   = clang++
  CXX11 = $(CXX) -std=c++11 
  OMP   =
endif

# ------------------------
# - Build specification. -
# ------------------------

# You shouldn't normally need to edit anything below here.

CXXFLAGS       = -W -Wno-unused-parameter -g -O3
CXXFLAGS_DEBUG = -g3 -fno-inline -O0 -Wall -Wextra 
BOOST_INC = -Iboost
BOOST_LIB = boost/stage/lib/libboost_program_options.a boost/stage/lib/libboost_random.a
LEMON_INC = -Ilemon/build -Ilemon/lemon-main-473c71baff72
LEMON_LIB = lemon/build/lemon/libemon.a -lpthread
CITY_INC  = -Icity/install/include
CITY_LIB  = city/install/lib/libcityhash.a
SH_INC 	  = -Isparsehash/include
DL_INC 	  = -Ideweylab
INC = $(BOOST_INC) $(LEMON_INC) $(CITY_INC) $(SH_INC) $(DL_INC)
LIB = $(BOOST_LIB) $(LEMON_LIB) $(CITY_LIB)
TEST_LIB = boost/stage/lib/libboost_unit_test_framework.a

.PHONY: all
all: ref-eval

.PHONY: debug
debug: CXXFLAGS += $(CXXFLAGS_DEBUG)
debug: ref-eval

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
	$(CXX) $(OMP) $(CXXFLAGS) $(INC) ref-eval.cpp $(LIB) -o ref-eval

all_tests := test_lazycsv test_line_stream test_blast test_psl test_pairset test_mask test_alignment_segment test_re_matched

.PHONY: test
test: test_msg ${all_tests} boost/finished lemon/finished city/finished sparsehash/finished
	./test_lazycsv
	./test_line_stream
	./test_blast
	./test_psl
	./test_pairset --show_progress
	./test_mask
	./test_alignment_segment
	./test_re_matched

.PHONY: test_msg
test_msg:
	@echo 
	@echo -------------------------------------------
	@echo - Building and/or running REF-EVAL tests. -
	@echo -------------------------------------------
	@echo 

test_lazycsv: test_lazycsv.cpp
	$(CXX) $(CXXFLAGS) $(INC) test_lazycsv.cpp $(LIB) $(TEST_LIB) -o test_lazycsv

test_line_stream: test_line_stream.cpp
	$(CXX) $(CXXFLAGS) $(INC) test_line_stream.cpp $(LIB) $(TEST_LIB) -o test_line_stream

test_blast: test_blast.cpp
	$(CXX) $(CXXFLAGS) $(INC) test_blast.cpp $(LIB) $(TEST_LIB) -o test_blast

test_psl: test_psl.cpp
	$(CXX) $(CXXFLAGS) $(INC) test_psl.cpp $(LIB) $(TEST_LIB) -o test_psl

test_pairset: test_pairset.cpp
	$(CXX) $(CXXFLAGS) $(INC) test_pairset.cpp $(LIB) $(TEST_LIB) -o test_pairset

test_mask: test_mask.cpp
	$(CXX) $(CXXFLAGS) $(INC) test_mask.cpp $(LIB) $(TEST_LIB) -o test_mask

test_alignment_segment: test_alignment_segment.cpp alignment_segment.hh
	$(CXX11) $(CXXFLAGS) $(INC) test_alignment_segment.cpp $(LIB) $(TEST_LIB) -o test_alignment_segment

test_re_matched: test_re_matched.cpp re_matched.hh
	$(CXX11) $(CXXFLAGS) $(INC) test_re_matched.cpp $(LIB) $(TEST_LIB) -o test_re_matched

.PHONY: clean
top-clean:
	-rm -f ref-eval ${all_tests}
	-rm -rf *.dSYM

.PHONY: clean
clean: top-clean
	-cd boost && $(MAKE) clean
	-cd lemon && $(MAKE) clean
	-cd city && $(MAKE) clean
	-cd sparsehash && $(MAKE) clean
