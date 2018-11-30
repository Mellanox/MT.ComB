/*
 * Copyright (c) 2017-2018 Mellanox Technologies Ltd. ALL RIGHTS RESERVED.
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "generic.h"

int is_mpi_based_test()
{
    return 1;
}

void special_usage(char *cmd) {
    fprintf(stderr, "\t-Dthrds\tDisable threaded support (call MPI_Init) (default: %s)\n",
            (want_thr_support)  ? "MPI_Init_thread" : "MPI_Init");
    fprintf(stderr, "\t-d\tUse separate communicator for each thread (default: %s)\n",
            (dup_comm) ? "enabled" : "disabled");
}

int special_process_args(char *optarg) {
    return 0;
}

int init_ctx() {
}

int connect_eps() {
}

int mem_map() { /* nothing to do */ }

void setup_thread_info_single(struct thread_info *ti)
{
    ti[0].tid = 0;
    if (dup_comm) {
        MPI_Comm_dup(MPI_COMM_WORLD, &ti[0].comm);
    } else {
        ti[0].comm = MPI_COMM_WORLD;
    }
}

void setup_thread_info_multi(struct thread_info *ti, int i)
{
    ti[i].tid = i;
    if (dup_comm) {
        MPI_Comm_dup(MPI_COMM_WORLD, &ti[i].comm);
    } else {
        ti[i].comm = MPI_COMM_WORLD;
    }
}

void cleanup_thread_info(struct thread_info *ti, int size)
{
    int i;
    if (dup_comm) {
        for (i = 0; i < size; i++) {
            MPI_Comm_free(&ti[i].comm);
        }
    }
}

void cleanup_ctx() {
}

void *init_sreqs(int num_sreqs, void *data_buf, int msg_sz, int tag) {
    MPI_Request *sreqs = calloc(num_sreqs, sizeof(MPI_Request));
    return (void *)sreqs;
}

void *init_rreqs(int num_rreqs, void *data_buf, int msg_sz, int tag) {
    MPI_Request *rreqs = calloc(num_rreqs, sizeof(MPI_Request));
    return (void *)rreqs;
}

void cleanup_sreqs(void *sreqs) {
    free(sreqs);
}

void cleanup_rreqs(void *rreqs) {
    free(rreqs);
}

void nonblocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *sreqs, int idx) {
    MPI_Isend(data_buf, msg_sz, MPI_BYTE, my_partner, tag, tinfo->comm, &((MPI_Request *)sreqs)[idx]);
}

void nonblocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *rreqs, int idx) {
    MPI_Irecv(data_buf, msg_sz, MPI_BYTE, my_partner, tag, tinfo->comm, &((MPI_Request *)rreqs)[idx]);
}

void wait_all_sreqs(void *sreqs, struct thread_info *tinfo, int num_sreqs) {
    MPI_Waitall(num_sreqs, (MPI_Request *)sreqs, MPI_STATUS_IGNORE);
}

void wait_all_rreqs(void *rreqs, struct thread_info *tinfo, int num_rreqs) {
    MPI_Waitall(num_rreqs, (MPI_Request *)rreqs, MPI_STATUS_IGNORE);
}

void blocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag) {
    MPI_Send(data_buf, msg_sz, MPI_BYTE, my_partner, tag, tinfo->comm);
}

void blocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag) {
    MPI_Recv(data_buf, msg_sz, MPI_BYTE, my_partner, tag, tinfo->comm, MPI_STATUS_IGNORE);
}

void wireup_progress(struct thread_info *tinfo){}
