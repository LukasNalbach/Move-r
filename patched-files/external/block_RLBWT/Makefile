
ifndef RUN_COUNT
RUN_COUNT = 32
endif

ifndef LARGE_BLOCK_SIZE
LARGE_BLOCK_SIZE = 16384
endif

ifndef SMALL_BLOCK_SIZE
SMALL_BLOCK_SIZE = 4096
endif

CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

CFLAGS = -std=c++2a -Wall -Wextra -Wshadow -pedantic -march=native -DLARGE_BLOCK_SIZE=$(LARGE_BLOCK_SIZE) \
         -DSMALL_BLOCK_SIZE=$(SMALL_BLOCK_SIZE) -DRUN_COUNT=$(RUN_COUNT) -DCACHE_LINE=$(CL)

HEADERS = include/reader.hpp include/block_rlbwt.hpp include/b_heap.hpp\
          include/byte_block.hpp include/byte_alphabet.hpp include/super_block.hpp \
		  include/types.hpp include/two_byte_block.hpp include/custom_alphabet.hpp \
		  include/one_byte_block.hpp include/d_block.hpp include/acgtn_alphabet.hpp \
		  include/alphabet.hpp include/vbyte_runs.hpp include/run_rlbwt.hpp 

.PHONY: clean update_git debug all

.DEFAULT: all

%/%.hpp:

all: gpp make_alphabet_header

gpp: make_bwt bench_bwt count_matches

make_bwt: make_bwt.cpp $(HEADERS)
	g++ $(CFLAGS) -DNDEBUG -Ofast -fpermissive -o make_bwt make_bwt.cpp -ldl

bench_bwt: bench_bwt.cpp $(HEADERS)
	g++ $(CFLAGS) -DNDEBUG -Ofast -fpermissive -o bench_bwt bench_bwt.cpp

make_alphabet_header: make_alphabet_header.cpp include/reader.hpp
	g++ $(CFLAGS) -DNDEBUG -Ofast -fpermissive -o make_alphabet_header make_alphabet_header.cpp

make_test_data: make_test_data.cpp $(HEADERS)
	g++ $(CFLAGS) -DNDEBUG -Ofast -fpermissive -o make_test_data make_test_data.cpp

count_matches: count_matches.cpp $(HEADERS)
	g++ $(CFLAGS) -DNDEBUG -Ofast -fpermissive -o count_matches count_matches.cpp

new_queries: new_queries.cpp $(HEADERS)
	g++ $(CFLAGS) -DNDEBUG -Ofast -o new_queries new_queries.cpp

debug: make_bwt.cpp bench_bwt.cpp $(HEADERS)
	g++ $(CFLAGS) -DDEBUG -g -o make_bwt make_bwt.cpp
	g++ $(CFLAGS) -DDEBUG -g -o bench_bwt bench_bwt.cpp
	g++ $(CFLAGS) -DDEBUG -g -o count_matches count_matches.cpp

clean:
	rm -f make_bwt bench_bwt count_matches make_alphabet_header count_matches make_test_data
