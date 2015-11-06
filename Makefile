#  This file is part of SSCA1.
#
#  Copyright (C) 2008-2015, UT-Battelle, LLC.
#
#  This product includes software produced by UT-Battelle, LLC under Contract No.
#  DE-AC05-00OR22725 with the Department of Energy.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the New BSD 3-clause software license (LICENSE).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  LICENSE for more details.
#
#  For more information please contact the SSCA1 developers at:
#  bakermb@ornl.gov


#valid compilers: gcc; pgcc, pathcc
#tested with gcc 4.3.3; pgi 8.0.6; pathscale 3.2
COMPILER=mpicc

#comment out this line to turn off debugging code
#DEBUG=1

ifeq ($(COMPILER), sgi)
        CC=cc
        COMMON_CFLAGS=-std=c99 -Wall -pipe -g -DUSE_SHMEM -DSGI_SHMEM
        OPTIMIZED_CFLAGS=-O3 -pipe -frename-registers
        DEBUG_CFLAGS=-O0 -ggdb -DDEBUG
	EXTRA_LIBS=-lsma -lmpi
endif

ifeq ($(COMPILER), cray)
        CC=cc
        COMMON_CFLAGS=-std=c99 -Wall -pipe -g -DUSE_SHMEM
        OPTIMIZED_CFLAGS=-O3 -pipe -frename-registers #-DUSE_PREFETCH
        DEBUG_CFLAGS=-O0 -ggdb -DDEBUG
endif

ifeq ($(COMPILER), pgi)
	CC=pgcc
	COMMON_CFLAGS=-c99
	OPTIMIZED_CFLAGS=-fastsse
	DEBUG_CFLAGS=-O0 -g
endif

ifeq ($(COMPILER), gcc)
	CC=gcc
	COMMON_CFLAGS=-std=c99 -Wall -pipe -g #-pg #-Werror #-pthread #-fopenmp
	OPTIMIZED_CFLAGS=-O3 -pipe -frename-registers -fopenmp
	DEBUG_CFLAGS=-O0 -ggdb -DDEBUG
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

ifeq ($(COMPILER), mpicc)
	CC=cc
	COMMON_CFLAGS=-std=c99 -Wall -pipe -g -DUSE_MPI3
	OPTIMIZED_CFLAGS=-O3 -pipe -frename-registers -DUSE_PREFETCH
	DEBUG_CFLAGS=-O0 -ggdb -DDEBUG
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

COMMON_OBJS=parameters.o gen_sim_matrix.o gen_scal_data.o pairwise_align.o glibc_sort.o scan_backwards.o util.o
MAIN_OBJS=main.o

all: $(MAIN_OBJS) $(COMMON_OBJS)
	$(LD) $(MAIN_OBJS) $(COMMON_OBJS) $(CFLAGS) $(LDFLAGS) -o ssca1 $(LIBS) $(EXTRA_LIBS)

clean:
	rm -f $(COMMON_OBJS) $(TEST_OBJS) $(MAIN_OBJS) *~ core.* ssca1 ssca1-test
