/*
 * Copyright (c) 2017-2018 Mellanox Technologies Ltd. ALL RIGHTS RESERVED.
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#ifndef GENERIC_H_
#define GENERIC_H_

#ifndef __USE_GNU
#define __USE_GNU
#endif

#define _GNU_SOURCE

#include <mpi.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sched.h>
#include <assert.h>

//#define DEBUG

#define WMB() \
{ asm volatile ("" : : : "memory"); }

#define MAX_HOSTS 2
#define MAX_HOSTHAME 1024

#define MODE_MT_SINGLE 0
#define MODE_MT_BASIC 1

struct thread_info {
    int tid;
    MPI_Comm comm;
    int worker_idx;
};

int is_mpi_based_test();

void *(*worker)(void *);
void *worker_b(void *);
void *worker_nb(void *);

void *init_sreqs(int num_sreqs, void *data_buf, int msg_sz, int tag);
void *init_rreqs(int num_rreqs, void *data_buf, int msg_sz, int tag);
void cleanup_sreqs(void *sreqs);
void cleanup_rreqs(void *rreqs);

void nonblocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *sreqs, int idx);
void nonblocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *rreqs, int idx);
void blocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag);
void blocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag);
void wait_all_sreqs(void *sreqs, struct thread_info *tinfo, int num_sreqs);
void wait_all_rreqs(void *rreqs, struct thread_info *tinfo, int num_rreqs);
void wireup_progress(struct thread_info *tinfo);
void print_results(MPI_Comm comm);
void allocate_global_buf();
void verify_buf();

int init_ctx();
void cleanup_ctx();
int connect_eps();

void setup_thread_info_single(struct thread_info *ti);
void setup_thread_info_multi(struct thread_info *ti, int idx);
void cleanup_thread_info(struct thread_info *ti, int size);
int mem_map();

void set_default_args();
void pre_scan_args(int argc, char **argv);
void special_usage(char *cmd);

void sync_worker(struct thread_info *tinfo);
void sync_master();
void usage(char *cmd);
int check_unsigned(char *str);
void process_args(int argc, char **argv);
int special_process_args(char *optarg);
void setup_binding(MPI_Comm comm);
void bind_worker(int tid);
MPI_Comm split_to_pairs();

enum bind_errcodes {
    BIND_OK,
    BIND_OVERSUBS,
    BIND_OVERLAP,
    BIND_GETAFF
};

extern int threads;
extern int win_size;
extern int iterations;
extern int warmup;
extern int dup_comm;
extern int msg_size;
extern int want_thr_support;
extern int intra_node;
extern int verify_mode;

extern int my_host_idx, my_rank_idx, my_leader;
extern int my_partner, i_am_sender;
extern int send_alignment;
extern int recv_alignment;
extern int buf_unit_size;
extern int remote_buf_unit_size;

extern char *global_buf;
extern double *results;

extern int volatile sync_cur_step;
extern int volatile *sync_thread_ready;
extern int volatile sync_start_all;

#endif /* GENERIC_H */
