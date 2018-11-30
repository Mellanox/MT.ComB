#
#  Copyright (c) 2017-2018 Mellanox Technologies Ltd. ALL RIGHTS RESERVED.
#
#   $COPYRIGHT$
#
#   Additional copyrights may follow
#
#   $HEADER$
#

CFLAGS ?= -g -lpthread
MPICC ?= 

UCX_HOME ?= 
UCX_CFLAGS = -I$(UCX_HOME)/include/
UCX_LDFLAGS = -L$(UCX_HOME)/lib/ -lucp -lucs -luct -lucm
MXM_HOME ?=
MXM_CFLAGS = -I$(MXM_HOME)/include/
MXM_LDFLAGS = -L$(MXM_HOME)/lib/ -lmxm 
LIBF_HOME ?=



all: mtcomb_mpip2p mtcomb_mpiosc mtcomb_ucxp2p mtcomb_ucxamo

mtcomb_mpip2p: generic.h generic.c mpi.c
	$(MPICC) $(CFLAGS) -o mtcomb_mpi generic.h generic.c mpi.c timeline.c

mtcomb_mpiosc: generic.h generic.c mpi_osc.c
	$(MPICC) $(CFLAGS) -o mtcomb_mpiosc generic.h generic.c mpi_osc.c

mtcomb_ucxp2p: generic.h generic.c ucx.c
	$(MPICC) -Wl,-rpath,$(UCX_HOME)/lib/ $(CFLAGS) $(UCX_CFLAGS) -o mtcomb_ucxp2p generic.h generic.c ucx.c $(UCX_LDFLAGS)

mtcomb_ucxamo: generic.h generic.c ucx.c
	$(MPICC) -Wl,-rpath,$(UCX_HOME)/lib/ $(CFLAGS) $(UCX_CFLAGS) -o mtcomb_ucxamo generic.h generic.c ucx_atomic.c $(UCX_LDFLAGS)

mtcomb_mxm: generic.h generic.c mxm.c
	$(MPICC) -Wl,-rpath,$(MXM_HOME)/lib/ $(CFLAGS) $(MXM_CFLAGS) -o mxm_bench generic.h generic.c mxm.c $(MXM_LDFLAGS)

#libf: generic.h generic.c libfabric.c
#	$(MPICC) -Wl,-rpath,$(LIBF_HOME)/lib/ -I$(LIBF_HOME)/include/ -L$(LIBF_HOME)/lib/ $(CFLAGS) -o libf_bench generic.h generic.c libfabric.c

clean:
	rm -f mtcomb_* *.o *.a
