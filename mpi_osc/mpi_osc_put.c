/*
 * Copyright (c) 2017-2018 Mellanox Technologies Ltd. ALL RIGHTS RESERVED.
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "../generic.h"

MPI_Win win;

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
}

int init_ctx() {
    /* TODO: add support for multi-comm case (multi-win)
     * in this case need to move it to the setup_thread* functions
     */
    MPI_Win_create(global_buf, threads*msg_size, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_lock_all(0, win);

}

int connect_eps() {
}

int mem_map() { /* nothing to do */ }

/* Same as for MPI */
void setup_thread_info_single(struct thread_info *ti)
{
    if (dup_comm) {
        MPI_Comm_dup(MPI_COMM_WORLD, &ti[0].comm);
    } else {
        ti[0].comm = MPI_COMM_WORLD;
    }
}

void setup_thread_info_multi(struct thread_info *ti, int i)
{
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
    MPI_Win_unlock_all(win);
    MPI_Win_free(&win);
}

void *init_sreqs(int num_sreqs, void *data_buf, int msg_sz, int tag) {
    return NULL;
}

void *init_rreqs(int num_rreqs, void *data_buf, int msg_sz, int tag) {
    return NULL;
}

void cleanup_sreqs(void *sreqs) {
}

void cleanup_rreqs(void *rreqs) {
}

void nonblocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *sreqs, int idx) {
    int delay = 0;
    while (delay) {
        sleep(1);
    }
    MPI_Put(data_buf, msg_sz, MPI_BYTE, my_partner, tinfo->tid * remote_buf_unit_size, msg_sz, MPI_BYTE, win);
}

void nonblocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *rreqs, int idx) {
}

void wait_all_sreqs(void *sreqs, struct thread_info *tinfo, int num_sreqs) {
    MPI_Win_flush(my_partner, win);
    // MPI_Win_flush(win);
}

void wait_all_rreqs(void *rreqs, struct thread_info *tinfo, int num_rreqs) {
}

void wireup_progress(struct thread_info *tinfo) {
}

void blocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag){
    MPI_Put(data_buf, msg_sz, MPI_BYTE, my_partner, tinfo->tid * remote_buf_unit_size, msg_sz, MPI_BYTE, win);
    MPI_Win_flush(my_partner, win);
}

void blocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag){
}
