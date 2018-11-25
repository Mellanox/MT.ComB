
CFLAGS ?= -g -lpthread
MPICC ?= /labhome/xinz/workplace/ompi/source/0debug/install/bin/mpicc

UCX_HOME ?= /labhome/xinz/workplace/ucx/install
MXM_HOME ?=
LIBF_HOME ?=

all: mpi mpi_osc ucx mxm libf ucx_atomic

mpi: generic.h generic.c mpi.c
	$(MPICC) $(CFLAGS) -o mpi_bench generic.h generic.c mpi.c

mpi_osc: generic.h generic.c mpi_osc.c
	$(MPICC) $(CFLAGS) -o mpi_osc_bench generic.h generic.c mpi_osc.c

ucx: generic.h generic.c ucx.c
	$(MPICC) -Wl,-rpath,$(UCX_HOME)/lib/ -I$(UCX_HOME)/include/ -L$(UCX_HOME)/lib/ $(CFLAGS) -lucp -lucs -o ucx_bench generic.h generic.c ucx.c

ucx_atomic: generic.h generic.c ucx.c
	$(MPICC) -Wl,-rpath,$(UCX_HOME)/lib/ -I$(UCX_HOME)/include/ -L$(UCX_HOME)/lib/ $(CFLAGS) -lucp -lucs -o ucx_atomic_bench generic.h generic.c ucx_atomic.c

mxm: generic.h generic.c mxm.c
	$(MPICC) -Wl,-rpath,$(MXM_HOME)/lib/ -I$(MXM_HOME)/include/ -L$(MXM_HOME)/lib/ $(CFLAGS) -lmxm -o mxm_bench generic.h generic.c mxm.c

libf: generic.h generic.c libfabric.c
	$(MPICC) -Wl,-rpath,$(LIBF_HOME)/lib/ -I$(LIBF_HOME)/include/ -L$(LIBF_HOME)/lib/ $(CFLAGS) -o libf_bench generic.h generic.c libfabric.c

clean:
	rm -f *_bench *.o *.a
