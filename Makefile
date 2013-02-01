CC = /usr/bin/g++
#DEBUG = -g3 -fno-inline -O0
#CFLAGS = -O3 -W -Wall -Wextra $(DEBUG)
DEBUG =
CFLAGS = -O3 -fopenmp -W -Wall -Wextra $(DEBUG)
LFLAGS = -Wall $(DEBUG)
INCLUDE = -I/Users/nathanae/Documents/fall12/cd/libdeweylab-git-install/include
LIBS = -lboost_program_options-mt -lboost_random-mt -L/Users/nathanae/Documents/fall12/cd/libdeweylab-git-install/lib -ldeweylab
TEST_LIBS = /opt/local/lib/libboost_unit_test_framework-mt.dylib

summarize: summarize.cpp
	$(CC) $(CFLAGS) $(INCLUDE) summarize.cpp $(LIBS) -o summarize

test: test_lazycsv test_line_stream test_blast test_pslx test_pairset test_mask
	./test_lazycsv
	./test_line_stream
	./test_blast
	./test_pslx
	./test_pairset
	./test_mask

test_lazycsv: test_lazycsv.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_lazycsv.cpp $(LIBS) $(TEST_LIBS) -o test_lazycsv

test_line_stream: test_line_stream.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_line_stream.cpp $(LIBS) $(TEST_LIBS) -o test_line_stream

test_blast: test_blast.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_blast.cpp $(LIBS) $(TEST_LIBS) -o test_blast

test_pslx: test_pslx.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_pslx.cpp $(LIBS) $(TEST_LIBS) -o test_pslx

test_pairset: test_pairset.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_pairset.cpp $(LIBS) $(TEST_LIBS) -o test_pairset

test_mask: test_mask.cpp
	$(CC) $(CFLAGS) $(INCLUDE) test_mask.cpp $(LIBS) $(TEST_LIBS) -o test_mask

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
