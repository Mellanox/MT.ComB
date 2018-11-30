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

void set_default_args() {
    worker = worker_nb;
    threads = 1;
    win_size = 256;
    iterations = 100;
    warmup = 10;
    dup_comm = 0;
    msg_size = 0;
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

int main(int argc, char *argv[]) {
    int rank, size, i, mt_level_act;
    char *sthreaded_env;
    pthread_t *id;
    MPI_Comm comm;

    set_default_args();

    /* unfortunately this is hackish */
    pre_scan_args(argc, argv);

    if (want_thr_support) {
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mt_level_act);
    } else {
        MPI_Init(&argc, &argv);
        mt_level_act = MPI_THREAD_SINGLE;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (want_thr_support && mt_level_act != MPI_THREAD_MULTIPLE) {
        if (rank == 0) {
            fprintf(stderr, "NOTE: no thread support!\n");
        }
    }

    process_args(argc, argv);

    if (threads > 1 && mt_level_act != MPI_THREAD_MULTIPLE) {
        if (rank == 0) {
            fprintf(stderr, "ERROR: %d threads requested but MPI implementation doesn't support THREAD_MULTIPLE\n", threads);
        }
        MPI_Finalize();
        exit(1);
    }

    comm = split_to_pairs();

    struct thread_info *ti = calloc(threads, sizeof(struct thread_info));
    results = calloc(threads, sizeof(double));
    id = calloc(threads, sizeof(*id));
    sync_thread_ready = calloc(threads, sizeof(int));

    /* allocate data buf */
    allocate_global_buf();

    MPI_Barrier(MPI_COMM_WORLD);

    if (threads == 1) {
        ti[0].tid = 0;
        if (dup_comm) {
            MPI_Comm_dup(MPI_COMM_WORLD, &ti[0].comm);
        } else {
            ti[0].comm = MPI_COMM_WORLD;
        }
        sync_start_all = sync_cur_step;

        worker((void*)ti);
    } else {
        /* Create the zero'ed array of ready flags for each thread */
        WMB();

        /* setup and create threads */
        for (i = 0; i < threads; i++) {
            ti[i].tid = i;
            if (dup_comm) {
                MPI_Comm_dup(MPI_COMM_WORLD, &ti[i].comm);
            } else {
                ti[i].comm = MPI_COMM_WORLD;
            }
            pthread_create(&id[i], NULL, worker, (void *) &ti[i]);
        }

        sync_master();

        /* wait for the test to finish */
        for (i = 0; i < threads; i++)
            pthread_join(id[i], NULL);
    }

    print_results(comm);

    if (verify_mode)
        verify_buf();

    if (dup_comm) {
        for (i = 0; i < threads; i++) {
            MPI_Comm_free(&ti[i].comm);
        }
    }

    free(id);
    free(results);
    free(global_buf);
    free(ti);

    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
