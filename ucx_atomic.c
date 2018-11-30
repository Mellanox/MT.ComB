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
#define UCX_MODE_MT_CTX_SHARED_RECV_WKR_SHARED 4

ucp_ep_h *eps = NULL;
ucp_worker_h *workers = NULL;
ucp_context_h *contexts = NULL;
ucp_rkey_h *rkeys = NULL;
ucp_mem_h *memhs = NULL;
void *remote_addr = NULL;

int num_workers = 1;
int num_contexts = 1;
int num_eps = 1;
int num_rkeys = 1;

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
    fprintf(stderr, "\t-Esctx_rswkr\tUse context-level threaded support, and on remote side use worker-level threaded support in UCX\n");
    fprintf(stderr, "\t-Eswkr\tUse worker-level threaded support in UCX\n");
}

int special_process_args(char *optarg) {
    if (!strcmp(optarg, "mctx")) {
        want_thr_support = UCX_MODE_MT_CTX_SEPARATED;
    } else if (!strcmp(optarg, "sctx")) {
        want_thr_support = UCX_MODE_MT_CTX_SHARED;
    } else if (!strcmp(optarg, "swkr")) {
        want_thr_support = UCX_MODE_MT_WKR_SHARED;
    } else if (!strcmp(optarg, "sctx_rswkr")) {
        want_thr_support = UCX_MODE_MT_CTX_SHARED_RECV_WKR_SHARED;
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
        while(req->completed == 0) {
            ucp_worker_progress(workers[worker_idx]);
        }
        request_release(req);
    }
}

void empty_function(void *request, ucs_status_t status) {
    struct request *req = (struct request *)request;
    assert(req->completed == 0);
    req->completed = 1;
}

ucs_status_t blocking_ep_flush(ucp_ep_h ep, ucp_worker_h worker)
{
    void *request = ucp_ep_flush_nb(ep, 0, empty_function);
    if (request == NULL) {
        return UCS_OK;
    } else if (UCS_PTR_IS_ERR(request)) {
        return UCS_PTR_STATUS(request);
    } else {
        ucs_status_t status;
        do {
            ucp_worker_progress(worker);
            status = ucp_request_check_status(request);
        } while (status == UCS_INPROGRESS);
        assert(((struct request *)request)->completed);
        request_release(request);
        return status;
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
    } else if (want_thr_support == UCX_MODE_MT_CTX_SHARED_RECV_WKR_SHARED) {
        if (i_am_sender) {
            params.mt_workers_shared = 1;
            worker_params.thread_mode = UCS_THREAD_MODE_SINGLE;
            num_workers = threads;
        } else {
            params.mt_workers_shared = 0;
            worker_params.thread_mode = UCS_THREAD_MODE_MULTI;
        }
    } else { /* no thread support */
        params.mt_workers_shared = 0;
        worker_params.thread_mode = UCS_THREAD_MODE_SINGLE;
    }

    params.field_mask |= UCP_PARAM_FIELD_FEATURES |
                         UCP_PARAM_FIELD_REQUEST_INIT |
                         UCP_PARAM_FIELD_REQUEST_SIZE;
    params.features = UCP_FEATURE_TAG | UCP_FEATURE_RMA;
    params.request_init = request_init;
    params.request_size = sizeof(struct request);

    /* create UCX contexts */
    contexts = calloc(num_contexts, sizeof(ucp_context_h));
    for (i = 0; i < num_contexts; i++) {
        ucp_config_t *config;
        ucp_config_read(NULL, NULL, &config);
        status = ucp_init(&params, config, &contexts[i]);
        assert(status == UCS_OK);
        if (want_thr_support == UCX_MODE_MT_CTX_SHARED ||
            (i_am_sender && want_thr_support == UCX_MODE_MT_CTX_SHARED_RECV_WKR_SHARED)) {
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
        if (want_thr_support == UCX_MODE_MT_WKR_SHARED ||
            (!i_am_sender && want_thr_support == UCX_MODE_MT_CTX_SHARED_RECV_WKR_SHARED)) {
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
    ucp_address_t *my_addr, *pair_addr;
    size_t my_addr_len, pair_addr_len;
    MPI_Status mpi_status;

    if (want_thr_support == UCX_MODE_MT_CTX_SHARED_RECV_WKR_SHARED) {
        if (!i_am_sender) {
            status = ucp_worker_get_address(workers[0], &my_addr, &my_addr_len);
            assert(status == UCS_OK);

            MPI_Send(&my_addr_len, sizeof(size_t), MPI_BYTE,
                     my_partner, 0, MPI_COMM_WORLD);
            MPI_Send(my_addr, my_addr_len, MPI_BYTE,
                     my_partner, 0, MPI_COMM_WORLD);

            ucp_worker_release_address(workers[0], my_addr);

        } else {
            MPI_Recv(&pair_addr_len, sizeof(size_t), MPI_BYTE,
                     my_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            pair_addr = malloc(pair_addr_len);

            MPI_Recv(pair_addr, pair_addr_len, MPI_BYTE,
                     my_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            num_eps = num_workers;
            eps = calloc(num_eps, sizeof(ucp_ep_h));
            for (i = 0; i < num_eps; i++) {
                ucp_ep_params_t ep_params;
                ep_params.field_mask = UCP_EP_PARAM_FIELD_REMOTE_ADDRESS;
                ep_params.address = pair_addr;
                status = ucp_ep_create(workers[i], &ep_params, &eps[i]);
                if (status == UCS_ERR_UNREACHABLE) {
                    return 1;
                }
                assert(status == UCS_OK);
            }

            free(pair_addr);
        }
    } else {
        num_eps = num_workers;
        eps = calloc(num_eps, sizeof(ucp_ep_h));
        for (i = 0; i < num_eps; i++) {
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
    }

    return 0;
}

int mem_map() {
    ucp_mem_map_params_t mem_params;
    ucs_status_t status;
    int i;

    if (!i_am_sender) {
        memhs = calloc(num_contexts, sizeof(ucp_mem_h));
    }

    if (want_thr_support != UCX_MODE_MT_CTX_SEPARATED) {
        size_t rkey_buffer_size, pair_rkey_size;
        void *rkey_buffer = NULL, *pair_rkey_buffer = NULL;
        void *pair_addr = NULL;

        if (!i_am_sender) {
            memset(&mem_params, 0, sizeof(ucp_mem_map_params_t));
            mem_params.field_mask = UCP_MEM_MAP_PARAM_FIELD_ADDRESS |
                                    UCP_MEM_MAP_PARAM_FIELD_LENGTH;
            mem_params.length = threads * buf_unit_size;
            mem_params.address = global_buf;

            status = ucp_mem_map(contexts[0], &mem_params, &memhs[0]);

            ucp_rkey_pack(contexts[0], memhs[0], &rkey_buffer,
                          &rkey_buffer_size);

            MPI_Send(&rkey_buffer_size, sizeof(size_t), MPI_BYTE,
                     my_partner, 1, MPI_COMM_WORLD);
            MPI_Send(rkey_buffer, rkey_buffer_size, MPI_BYTE,
                     my_partner, 1, MPI_COMM_WORLD);

            ucp_rkey_buffer_release(rkey_buffer);

            MPI_Send(&global_buf, sizeof(void *), MPI_BYTE,
                     my_partner, 1, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&pair_rkey_size, sizeof(size_t), MPI_BYTE,
                     my_partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            pair_rkey_buffer = malloc(pair_rkey_size);

            MPI_Recv(pair_rkey_buffer, pair_rkey_size, MPI_BYTE,
                     my_partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            num_rkeys = num_workers;
            rkeys = calloc(num_rkeys, sizeof(ucp_rkey_h));
            for (i = 0; i < num_rkeys; i++) {
                status = ucp_ep_rkey_unpack(eps[i], pair_rkey_buffer,
                                            &rkeys[i]);
            }

            MPI_Recv(&remote_addr, sizeof(void *), MPI_BYTE,
                     my_partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        if (i_am_sender) {
            MPI_Recv(&remote_addr, sizeof(void *), MPI_BYTE,
                     my_partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            num_rkeys = num_workers;
            rkeys = calloc(num_rkeys, sizeof(ucp_rkey_h));

        } else {
            MPI_Send(&global_buf, sizeof(void *), MPI_BYTE,
                     my_partner, 1, MPI_COMM_WORLD);
        }

        for (i = 0; i < num_contexts; i++) {
            size_t rkey_buffer_size, pair_rkey_size;
            void *rkey_buffer = NULL, *pair_rkey_buffer = NULL;

            if (!i_am_sender) {
                memset(&mem_params, 0, sizeof(ucp_mem_map_params_t));
                mem_params.field_mask = UCP_MEM_MAP_PARAM_FIELD_ADDRESS | UCP_MEM_MAP_PARAM_FIELD_LENGTH;
                mem_params.length = threads * buf_unit_size;
                mem_params.address = global_buf;

                status = ucp_mem_map(contexts[i], &mem_params, &memhs[i]);

                ucp_rkey_pack(contexts[i], memhs[i], &rkey_buffer,
                              &rkey_buffer_size);

                MPI_Send(&rkey_buffer_size, sizeof(size_t), MPI_BYTE,
                         my_partner, 1, MPI_COMM_WORLD);
                MPI_Send(rkey_buffer, rkey_buffer_size, MPI_BYTE,
                         my_partner, 1, MPI_COMM_WORLD);

                ucp_rkey_buffer_release(rkey_buffer);
            } else {
                MPI_Recv(&pair_rkey_size, sizeof(size_t), MPI_BYTE,
                         my_partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                pair_rkey_buffer = malloc(pair_rkey_size);

                MPI_Recv(pair_rkey_buffer, pair_rkey_size, MPI_BYTE,
                         my_partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                status = ucp_ep_rkey_unpack(eps[i], pair_rkey_buffer,
                                            &rkeys[i]);
            }
        }
    }

    return 0;
}

void cleanup_ctx() {
    int i;

    if (memhs != NULL) {
        for (i = 0; i < num_contexts; i++) {
            ucp_mem_unmap(contexts[i], memhs[i]);
        }
    }

    if (rkeys != NULL) {
        for (i = 0; i < num_rkeys; i++) {
            ucp_rkey_destroy(rkeys[i]);
        }
        free(rkeys);
    }
    if (eps != NULL) {
        for (i = 0; i < num_eps; i++) {
            ucp_ep_destroy(eps[i]);
        }
        free(eps);
    }

    for (i = 0; i < num_workers; i++) {
        ucp_worker_destroy(workers[i]);
    }
    for (i = 0; i < num_contexts; i++) {
        ucp_cleanup(contexts[i]);
    }
    free(workers);
    free(contexts);
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
    ucp_put_nbi(eps[tinfo->worker_idx], data_buf, msg_sz, (uint64_t)remote_addr + tinfo->tid * remote_buf_unit_size, rkeys[tinfo->worker_idx]);
}

void nonblocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *rreqs, int idx) {
    // nothing
}

void blocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag) {
    ucp_put_nbi(eps[tinfo->worker_idx], data_buf, msg_sz, (uint64_t)remote_addr + tinfo->tid * remote_buf_unit_size, rkeys[tinfo->worker_idx]);
    blocking_ep_flush(eps[tinfo->worker_idx], workers[tinfo->worker_idx]);
}

void blocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag) {
    // nothing
}

void wait_all_sreqs(void *sreqs, struct thread_info *tinfo, int num_sreqs) {
    blocking_ep_flush(eps[tinfo->worker_idx], workers[tinfo->worker_idx]);
}

void wait_all_rreqs(void *rreqs, struct thread_info *tinfo, int num_rreqs) {
}

void wireup_progress(struct thread_info *tinfo) {
    if (threads == 1) {
        char syncbuf[4];
        struct request *sync_req;
        if (i_am_sender) {
            sync_req = (struct request *)ucp_tag_send_nb(eps[tinfo->worker_idx],
                                                         syncbuf, 4,
                                                         ucp_dt_make_contig(1),
                                                         0x1332,
                                                         send_callback);
        } else {
            sync_req = (struct request *)ucp_tag_recv_nb(workers[tinfo->worker_idx],
                                                         syncbuf, 4,
                                                         ucp_dt_make_contig(1),
                                                         0x1332, 0xffff,
                                                         recv_callback);
        }
        wait_and_release(tinfo->worker_idx, sync_req);
    } else {
        ucp_worker_progress(workers[tinfo->worker_idx]);
    }
}

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

    allocate_global_buf();

    if (ret = mem_map()) {
        MPI_Finalize();
        exit(1);
    }

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
            if (want_thr_support == UCX_MODE_MT_CTX_SHARED_RECV_WKR_SHARED) {
                if (i_am_sender) {
                    ti[i].worker_idx = i;
                } else {
                    ti[i].worker_idx = 0;
                }
            } else {
                ti[i].worker_idx = (want_thr_support != UCX_MODE_MT_WKR_SHARED ? i : 0);
            }
            pthread_create(&id[i], NULL, worker, (void *)&ti[i]);
        }

        sync_master();

        /* wait for the test to finish */
        for (i = 0; i < threads; i++)
            pthread_join(id[i], NULL);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    print_results(comm);

    if (verify_mode)
        verify_buf();

    free(id);
    free(results);
    free(ti);

    cleanup_ctx();

    free(global_buf);

    MPI_Comm_free(&comm);

    MPI_Finalize();

    return 0;
}
