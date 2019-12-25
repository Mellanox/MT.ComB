#
#  Copyright (c) 2017-2018 Mellanox Technologies Ltd. ALL RIGHTS RESERVED.
#
#   $COPYRIGHT$
#
#   Additional copyrights may follow
#
#   $HEADER$
#

CFLAGS ?= -g  -lpthread -DMULTI_MESSAGES
MPICC ?= mpicc

UCX_HOME ?= 
UCX_CFLAGS = -I$(UCX_HOME)/include/
UCX_LDFLAGS = -L$(UCX_HOME)/lib/ -lucp -lucs -luct -lucm
MXM_HOME ?=
MXM_CFLAGS = -I$(MXM_HOME)/include/
MXM_LDFLAGS = -L$(MXM_HOME)/lib/ -lmxm 
LIBF_HOME ?=

GENERIC_DEPS = generic.h generic.c mpi.c timeline.c
GENERIC_CFILES = generic.c timeline.c

all: mtcomb_mpip2p mtcomb_mpiosc mtcomb_ucxp2p mtcomb_ucxamo

mtcomb_mpip2p: $(GENERIC_DEPS) mpi.c
	$(MPICC) $(CFLAGS) -o mtcomb_mpi mpi.c $(GENERIC_CFILES)

mtcomb_mpiosc: $(GENERIC_DEPS) mpi_osc.c
	$(MPICC) $(CFLAGS) -o mtcomb_mpiosc mpi_osc.c $(GENERIC_CFILES)

mtcomb_ucxp2p: $(GENERIC_DEPS) ucx.c
	$(MPICC) -Wl,-rpath,$(UCX_HOME)/lib/ $(CFLAGS) $(UCX_CFLAGS) -o mtcomb_ucxp2p ucx.c $(GENERIC_CFILES) $(UCX_LDFLAGS)

mtcomb_ucxamo: $(GENERIC_DEPS) ucx.c
	$(MPICC) -Wl,-rpath,$(UCX_HOME)/lib/ $(CFLAGS) $(UCX_CFLAGS) -o mtcomb_ucxamo ucx_atomic.c $(GENERIC_CFILES) $(UCX_LDFLAGS)

mtcomb_mxm: $(GENERIC_DEPS) mxm.c
	$(MPICC) -Wl,-rpath,$(MXM_HOME)/lib/ $(CFLAGS) $(MXM_CFLAGS) -o mtcomb_mxm mxm.c $(GENERIC_CFILES) $(MXM_LDFLAGS)

#libf: generic.h generic.c libfabric.c
#	$(MPICC) -Wl,-rpath,$(LIBF_HOME)/lib/ -I$(LIBF_HOME)/include/ -L$(LIBF_HOME)/lib/ $(CFLAGS) -o libf_bench generic.h generic.c libfabric.c

clean:
	rm -f mtcomb_* *.o *.a *.eps timeline *.gpl
