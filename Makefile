#valid compilers: gcc; pgcc, pathcc
#tested with gcc 4.3.3; pgi 8.0.6; pathscale 3.2
COMPILER=gcc

#comment out this line to turn off debugging code
DEBUG=0

ifeq ($(COMPILER), pgi)
	CC=pgcc
	COMMON_CFLAGS=-c99
	OPTIMIZED_CFLAGS=-fastsse
	DEBUG_CFLAGS=-O0 -g
endif

ifeq ($(COMPILER), gcc)
	CC=gcc 
	COMMON_CFLAGS=-std=c99 -Wall -pipe -g -pthread #-fopenmp
	OPTIMIZED_CFLAGS=-O3 -pipe -frename-registers
	DEBUG_CFLAGS=-O0 -ggdb
endif

ifeq ($(COMPILER), pathscale)
	CC=pathcc
	COMMON_CFLAGS=-fullwarn -std=c99 -Wall
	OPTIMIZED_CFLAGS=-g0 -O2
	DEBUG_CFLAGS=-O0 -g2
endif

ifeq ($(COMPILER), intel)
        CC=icc
        COMMON_CFLAGS=-std=c99
        OPTIMIZED_CFLAGS=-fast -xhost
        DEBUG_CFLAGS=-O0 -g2
endif

ifndef DEBUG
	DEBUG=0
endif

ifneq ($(DEBUG), 0)
	CFLAGS=$(COMMON_CFLAGS) $(DEBUG_CFLAGS)
else
	CFLAGS=$(COMMON_CFLAGS) $(OPTIMIZED_CFLAGS)
endif

LD=$(CC)

CPPFLAGS=-I.
LIBS=-lm

COMMON_OBJS=parameters.o gen_sim_matrix.o gen_scal_data.o pairwise_align.o glibc_sort.o scan_backwards.o locate_similar.o util.o global_align.o multiple_align.o
MAIN_OBJS=main.o
TEST_OBJS=test_data.o test_main.o test_locate_similar.o merge_test.o test_multiple_align.o

all: $(MAIN_OBJS) $(COMMON_OBJS)
	$(LD) $(MAIN_OBJS) $(COMMON_OBJS) $(CFLAGS) $(LDFLAGS) -o ssca1 $(LIBS)

# this will make a test to compare the output of this version against the expected output of the Matlab version
test: $(TEST_OBJS) $(COMMON_OBJS)
	$(LD) $(TEST_OBJS) $(COMMON_OBJS) $(CFLAGS) $(LDFLAGS) -o ssca1-test $(LIBS)

clean:
	rm -f $(COMMON_OBJS) $(TEST_OBJS) $(MAIN_OBJS) *~ core.* ssca1 ssca1-test
