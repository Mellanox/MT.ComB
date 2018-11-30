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
#include <ucp/api/ucp.h>

#define UCX_MODE_MT_CTX_SEPARATED MODE_MT_BASIC
#define UCX_MODE_MT_CTX_SHARED 2
#define UCX_MODE_MT_WKR_SHARED 3

ucp_ep_h *eps = NULL;
ucp_worker_h *workers = NULL;
ucp_context_h *contexts = NULL;
int num_workers = 1;
int num_contexts = 1;

void set_default_args() {
    worker = worker_nb;
    threads = 1;
    win_size = 256;
    iterations = 50;
    warmup = 10;
    msg_size = 0;
}

void special_usage(char *cmd) {
    fprintf(stderr, "\t-Dthrds\tDisable threaded support in UCX\n");
    fprintf(stderr, "\t-Emctx\tSeparate context per thread (default)\n");
    fprintf(stderr, "\t-Esctx\tUse context-level threaded support in UCX\n");
    fprintf(stderr, "\t-Eswkr\tUse worker-level threaded support in UCX\n");
}

int special_process_args(char *optarg) {
    if (!strcmp(optarg, "mctx")) {
        want_thr_support = UCX_MODE_MT_CTX_SEPARATED;
    } else if (!strcmp(optarg, "sctx")) {
        want_thr_support = UCX_MODE_MT_CTX_SHARED;
    } else if (!strcmp(optarg, "swkr")) {
        want_thr_support = UCX_MODE_MT_WKR_SHARED;
    } else {
        return 1;
    }

    return 0;
}

struct request {
    int completed;
};

void request_init(void *request) {
    struct request *req = (struct request *)request;
    req->completed = 0;
}

void request_release(struct request *req) {
    req->completed = 0;
    ucp_request_release(req);
}

void send_callback(void *request, ucs_status_t status) {
    struct request *req = (struct request *)request;
    assert(req->completed == 0);
    req->completed = 1;
}

void recv_callback(void *request, ucs_status_t status, ucp_tag_recv_info_t *info) {
    struct request *req = (struct request *)request;
    assert(req->completed == 0);
    req->completed = 1;
}

void wait_and_release(int worker_idx, struct request *req) {
    if (req != NULL) {
        while (req->completed == 0) {
            ucp_worker_progress(workers[worker_idx]);
        }
        request_release(req);
    }
}

void waitall_and_release(int req_num, int worker_idx, struct request **reqs) {
    int i;
    for (i = 0; i < req_num; i++) {
        wait_and_release(worker_idx, reqs[i]);
    }
}

int init_ctx() {
    ucp_params_t params;
    ucp_worker_params_t worker_params;
    ucs_status_t status;
    int i;

    /* init ucx environment */
    params.field_mask = UCP_PARAM_FIELD_MT_WORKERS_SHARED;
    worker_params.field_mask = UCP_WORKER_PARAM_FIELD_THREAD_MODE;
    if (want_thr_support == UCX_MODE_MT_CTX_SEPARATED) {
        /* separate context per thread */
        params.mt_workers_shared = 0;
        worker_params.thread_mode = UCS_THREAD_MODE_SINGLE;
        num_contexts = threads;
        num_workers = threads;
    } else if (want_thr_support == UCX_MODE_MT_WKR_SHARED) {
        /* worker-level thread support */
        params.mt_workers_shared = 0;
        worker_params.thread_mode = UCS_THREAD_MODE_MULTI;
    } else if (want_thr_support == UCX_MODE_MT_CTX_SHARED) {
        /* context-level thread support */
        params.mt_workers_shared = 1;
        worker_params.thread_mode = UCS_THREAD_MODE_SINGLE;
        num_workers = threads;
    } else { /* no thread support */
        params.mt_workers_shared = 0;
        worker_params.thread_mode = UCS_THREAD_MODE_SINGLE;
    }

    params.field_mask |= UCP_PARAM_FIELD_FEATURES |
                         UCP_PARAM_FIELD_REQUEST_INIT |
                         UCP_PARAM_FIELD_REQUEST_SIZE;
    params.features = UCP_FEATURE_TAG;
    params.request_init = request_init;
    params.request_size = sizeof(struct request);

    /* create UCX contexts */
    contexts = calloc(num_contexts, sizeof(ucp_context_h));
    for (i = 0; i < num_contexts; i++) {
        ucp_config_t *config;
        ucp_config_read(NULL, NULL, &config);
        status = ucp_init(&params, config, &contexts[i]);
        assert(status == UCS_OK);
        if (want_thr_support == UCX_MODE_MT_CTX_SHARED) {
            ucp_context_attr_t attrs;
            attrs.field_mask = UCP_ATTR_FIELD_THREAD_MODE;
            status = ucp_context_query(contexts[i], &attrs);
            assert(status == UCS_OK);
            if (attrs.thread_mode != UCS_THREAD_MODE_MULTI) {
                fprintf(stderr, "NOTE: no thread support in UCP context!\n");
                return 1;
            }
        }
    }

    /* create UCX workers */
    workers = calloc(num_workers, sizeof(ucp_worker_h));
    for (i = 0; i < num_workers; i++) {
        int ctx_idx = (want_thr_support == UCX_MODE_MT_CTX_SEPARATED) ? i : 0;
        status = ucp_worker_create(contexts[ctx_idx], &worker_params, &workers[i]);
        assert(status == UCS_OK);
        if (want_thr_support == UCX_MODE_MT_WKR_SHARED) {
            ucp_worker_attr_t worker_attrs;
            worker_attrs.field_mask = UCP_WORKER_ATTR_FIELD_THREAD_MODE;
            status = ucp_worker_query(workers[i], &worker_attrs);
            assert(status == UCS_OK);
            if (worker_attrs.thread_mode != UCS_THREAD_MODE_MULTI) {
                fprintf(stderr, "NOTE: no thread support in UCP worker!\n");
                return 1;
            }
        }
    }

    return 0;
}

int connect_eps() {
    int i;
    ucs_status_t status;

    eps = calloc(num_workers, sizeof(ucp_ep_h));
    for (i = 0; i < num_workers; i++) {
        ucp_address_t *my_addr, *pair_addr;
        size_t my_addr_len, pair_addr_len;
        MPI_Status mpi_status;

        status = ucp_worker_get_address(workers[i], &my_addr, &my_addr_len);
        assert(status == UCS_OK);

        MPI_Sendrecv(&my_addr_len, sizeof(size_t), MPI_BYTE, my_partner, 0,
                     &pair_addr_len, sizeof(size_t), MPI_BYTE, my_partner, 0,
                     MPI_COMM_WORLD, &mpi_status);

        pair_addr = malloc(pair_addr_len);

        MPI_Sendrecv(my_addr, my_addr_len, MPI_BYTE, my_partner, 1,
                     pair_addr, pair_addr_len, MPI_BYTE, my_partner, 1,
                     MPI_COMM_WORLD, &mpi_status);

        ucp_ep_params_t ep_params;
        ep_params.field_mask = UCP_EP_PARAM_FIELD_REMOTE_ADDRESS;
        ep_params.address = pair_addr;
        status = ucp_ep_create(workers[i], &ep_params, &eps[i]);
        if (status == UCS_ERR_UNREACHABLE) {
            return 1;
        }
        assert(status == UCS_OK);

        ucp_worker_release_address(workers[i], my_addr);
        free(pair_addr);
    }

    return 0;
}

void cleanup_ctx() {
    int i;
    for (i = 0; i < num_workers; i++) {
        ucp_ep_destroy(eps[i]);
        ucp_worker_destroy(workers[i]);
    }
    for (i = 0; i < num_contexts; i++) {
        ucp_cleanup(contexts[i]);
    }
    free(eps);
    free(workers);
    free(contexts);
}

void *init_sreqs(int num_sreqs, void *data_buf, int msg_sz, int tag) {
    struct request **sreqs = calloc(num_sreqs, sizeof(struct request *));
    memset(sreqs, 0, num_sreqs * sizeof(struct request *));
    return (void *)sreqs;
}

void *init_rreqs(int num_rreqs, void *data_buf, int msg_sz, int tag) {
    struct request **rreqs = calloc(num_rreqs, sizeof(struct request *));
    memset(rreqs, 0, num_rreqs * sizeof(struct request *));
    return (void *)rreqs;
}

void cleanup_sreqs(void *sreqs) {
    free(sreqs);
}

void cleanup_rreqs(void *rreqs) {
    free(rreqs);
}

void nonblocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *sreqs, int idx) {
    struct request **ucx_reqs = (struct request **)sreqs;
    ucx_reqs[idx] = (struct request *)ucp_tag_send_nb(eps[tinfo->worker_idx],
                                                      data_buf, msg_sz,
                                                      ucp_dt_make_contig(1),
                                                      0x1337+tag,
                                                      send_callback);
}

void nonblocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *rreqs, int idx) {
    struct request **ucx_reqs = (struct request **)rreqs;
    ucx_reqs[idx] = (struct request *)ucp_tag_recv_nb(workers[tinfo->worker_idx],
                                                      data_buf, msg_sz,
                                                      ucp_dt_make_contig(1),
                                                      0x1337+tag, 0xffff,
                                                      recv_callback);
}

void blocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag) {
    struct request *sync_req = (struct request *)ucp_tag_send_nb(eps[tinfo->worker_idx],
                                                                 data_buf, msg_sz,
                                                                 ucp_dt_make_contig(1),
                                                                 0x1337+tag,
                                                                 send_callback);
    wait_and_release(tinfo->worker_idx, sync_req);
}

void blocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag) {
    struct request *sync_req = (struct request *)ucp_tag_recv_nb(workers[tinfo->worker_idx],
                                                                 data_buf, msg_sz,
                                                                 ucp_dt_make_contig(1),
                                                                 0x1337+tag, 0xffff,
                                                                 recv_callback);
    wait_and_release(tinfo->worker_idx, sync_req);
}

void wait_all_sreqs(void *sreqs, struct thread_info *tinfo, int num_sreqs) {
    char syncbuf[4];
    waitall_and_release(num_sreqs, tinfo->worker_idx, (struct request **)sreqs);
    blocking_recv(syncbuf, 4, tinfo, 123423);
}

void wait_all_rreqs(void *rreqs, struct thread_info *tinfo, int num_rreqs) {
    char syncbuf[4];
    waitall_and_release(num_rreqs, tinfo->worker_idx, (struct request **)rreqs);
    blocking_send(syncbuf, 4, tinfo, 123423);
}

void wireup_progress(struct thread_info *tinfo){}

int main(int argc, char *argv[]) {
    int rank, size, i;
    pthread_t *id;
    MPI_Comm comm;
    int ret, mt_level_act;

    set_default_args();

    pre_scan_args(argc, argv);

    if (!want_thr_support) {
        MPI_Init(&argc, &argv);
    } else {
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mt_level_act);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    process_args(argc, argv);

    if (threads > 1) {
        /* TODO: print out error message */
        assert(want_thr_support);
    }

    comm = split_to_pairs();

    if (ret = init_ctx()) {
        MPI_Finalize();
        exit(1);
    }

    if (ret = connect_eps()) {
        MPI_Finalize();
        exit(1);
    }

    struct thread_info *ti = calloc(threads, sizeof(struct thread_info));
    results = calloc(threads, sizeof(double));
    id = calloc(threads, sizeof(*id));
    sync_thread_ready = calloc(threads, sizeof(int));

    /* allocate data buf */
    allocate_global_buf();

    MPI_Barrier(MPI_COMM_WORLD);

    if (threads == 1) {
        ti[0].tid = 0;
        ti[0].worker_idx = 0;
        sync_start_all = sync_cur_step;
        worker((void*)ti);
    } else {
        /* Create the zero'ed array of ready flags for each thread */
        WMB();

        /* setup and create threads */
        for (i = 0; i < threads; i++) {
            ti[i].tid = i;
            ti[i].worker_idx = (want_thr_support != UCX_MODE_MT_WKR_SHARED ? i : 0);
            pthread_create(&id[i], NULL, worker, (void *)&ti[i]);
        }

        sync_master();

        /* wait for the test to finish */
        for (i = 0; i < threads; i++)
            pthread_join(id[i], NULL);
    }

    print_results(comm);

    if (verify_mode)
        verify_buf();

    free(id);
    free(results);
    free(global_buf);
    free(ti);

    cleanup_ctx();

    MPI_Comm_free(&comm);

    MPI_Finalize();

    return 0;
}
