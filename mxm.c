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
#include <mxm/api/mxm_api.h>

mxm_mq_h mq;
mxm_conn_h conn;
mxm_ep_h ep;
mxm_h mxmh;

int *quit_flag = NULL;

#define MXM_TEST_HID_QUIT 4

int is_mpi_based_test()
{
    return 0;
}

void special_usage(char *cmd) {
    fprintf(stderr, "\t-Dthrds\tDisable threaded support in MXM\n");
}

int special_process_args(char *optarg) {
    return 0;
}

static void init_req(mxm_req_base_t *req, void *data_addr, size_t data_len) {
    req->state = MXM_REQ_NEW;
    req->mq = mq;
    req->conn = conn;
    req->completed_cb = NULL;
    req->data_type = MXM_REQ_DATA_BUFFER;
    req->error = MXM_OK;
    req->data.buffer.ptr = data_addr;
    req->data.buffer.length = data_len;
}

static void init_send_req(mxm_send_req_t *sreq, void *data_addr, size_t data_len, int tag) {
    init_req(&sreq->base, data_addr, data_len);
    sreq->flags = 0;
    sreq->opcode = MXM_REQ_OP_SEND;
    sreq->op.send.tag = tag;
}

static void init_recv_req(mxm_recv_req_t *rreq, void *data_addr, size_t data_len, int tag) {
    init_req(&rreq->base, data_addr, data_len);
    rreq->tag = tag;
    rreq->tag_mask = -1;
}

int init_ctx() {
    mxm_context_opts_t *mxm_opts;
    mxm_ep_opts_t *ep_opts;
    mxm_error_t error;

    setenv("MXM_TLS", "ud,dc,rc", 0);
    setenv("MXM_IB_RX_QUEUE_LEN", "8192", 0);
    setenv("MXM_SINGLE_THREAD", want_thr_support == 0 ? "y" : "n", 1);

    error = mxm_config_read_opts(&mxm_opts, &ep_opts, NULL, NULL, 0);
    if (error != MXM_OK) {
        return 1;
    }

    error = mxm_init(mxm_opts, &mxmh);
    if (error != MXM_OK) {
        return 1;
    }

    error = mxm_ep_create_internal(mxmh, ep_opts, 1, &ep);
    if (error != MXM_OK) {
        return 1;
    }

    mxm_config_free_context_opts(mxm_opts);
    mxm_config_free_ep_opts(ep_opts);

    return 0;
}

int connect_eps() {
    char local_addr[256], remote_addr[256];
    size_t addr_len = 256;
    mxm_error_t error;
    MPI_Status mpi_status;

    memset(local_addr, 0, 256);
    memset(remote_addr, 0, 256);

    error = mxm_ep_get_address(ep, local_addr, &addr_len);
    if (error != MXM_OK) {
        return 1;
    }

    MPI_Sendrecv(local_addr, 256, MPI_BYTE, my_partner, 1,
                 remote_addr, 256, MPI_BYTE, my_partner, 1,
                 MPI_COMM_WORLD, &mpi_status);

    error = mxm_ep_connect(ep, remote_addr, &conn);
    if (error != MXM_OK) {
        return 1;
    }

    error = mxm_mq_create(mxmh, 0x1337, &mq);
    if (error != MXM_OK) {
        return 1;
    }

    error = mxm_ep_wireup(ep);
    if (error != MXM_OK) {
        return 1;
    }

    return 0;
}

int mem_map() { /* nothing to do */ }

void cleanup_ctx() {
    mxm_mq_destroy(mq);
    mxm_ep_disconnect(conn);
    mxm_ep_destroy(ep);
    mxm_cleanup(mxmh);
}

static void am_handler_quit(mxm_conn_h conn, mxm_imm_t imm, void *data, size_t len, size_t offset, int is_lf)
{
    assert(quit_flag != NULL);
    quit_flag[imm] = 1;
}

static mxm_error_t send_quit_am(int tid)
{
    mxm_send_req_t sreq;
    sreq.base.mq                  = mq;
    sreq.base.completed_cb        = NULL;
    sreq.base.data_type           = MXM_REQ_DATA_BUFFER;
    sreq.base.data.buffer.ptr     = NULL;
    sreq.base.data.buffer.length  = 0;
    sreq.opcode                   = MXM_REQ_OP_AM;
    sreq.flags                    = 0;
    sreq.op.am.hid                = MXM_TEST_HID_QUIT;
    sreq.op.am.imm_data           = tid;
    sreq.base.state               = MXM_REQ_NEW;
    sreq.base.conn                = conn;
    mxm_req_send(&sreq);
    mxm_req_wait(&sreq.base);
    return sreq.base.error;
}

void *init_sreqs(int num_sreqs, void *data_buf, int msg_sz, int tag) {
    int i;
    mxm_send_req_t *sreqs = calloc(num_sreqs, sizeof(mxm_send_req_t));
    for (i = 0; i < num_sreqs; i++) {
        init_send_req(&sreqs[i], (void *)data_buf, msg_sz, tag);
    }
    return (void *)sreqs;
}

void *init_rreqs(int num_rreqs, void *data_buf, int msg_sz, int tag) {
    int i;
    mxm_recv_req_t *rreqs = calloc(num_rreqs, sizeof(mxm_recv_req_t));
    for (i = 0; i < num_rreqs; i++) {
        init_recv_req(&rreqs[i], (void *)data_buf, msg_sz, tag);
    }
    return (void *)rreqs;
}

void cleanup_sreqs(void *sreqs) {
    free(sreqs);
}

void cleanup_rreqs(void *rreqs) {
    free(rreqs);
}

void nonblocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *sreqs, int idx) {
    mxm_req_send(&((mxm_send_req_t *)sreqs)[idx]);
}

void nonblocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *rreqs, int idx) {
    mxm_req_recv(&((mxm_recv_req_t *)rreqs)[idx]);
}

void blocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag) {
    mxm_send_req_t sreq;
    init_send_req(&sreq, (void *)data_buf, msg_sz, tag);
    mxm_req_send(&sreq);
    mxm_req_wait(&sreq.base);
}

void blocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag) {
    mxm_recv_req_t rreq;
    init_recv_req(&rreq, (void *)data_buf, msg_sz, tag);
    mxm_req_recv(&rreq);
    mxm_req_wait(&rreq.base);
}

void wait_all_sreqs(void *sreqs, struct thread_info *tinfo, int num_sreqs) {
    int i;
    for (i = 0; i < num_sreqs; i++) {
        mxm_req_wait(&((mxm_send_req_t *)sreqs)[i].base);
    }
}

void wait_all_rreqs(void *rreqs, struct thread_info *tinfo, int num_rreqs) {
    int i;
    for (i = 0; i < num_rreqs; i++) {
        mxm_req_wait(&((mxm_recv_req_t *)rreqs)[i].base);
    }
}

void wireup_progress(struct thread_info *tinfo){}

int main(int argc, char *argv[]) {
    int rank, size, i;
    pthread_t *id;
    MPI_Comm comm;
    int ret;

    set_default_args();

    pre_scan_args(argc, argv);

    MPI_Init(&argc, &argv);

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

    quit_flag = calloc(threads, sizeof(int));
    memset(quit_flag, 0, threads);
    mxm_set_am_handler(mxmh, MXM_TEST_HID_QUIT, am_handler_quit, 0);

    if (threads == 1) {
        ti[0].tid = 0;
        sync_start_all = sync_cur_step;
        worker((void*)ti);
    } else {
        /* Create the zero'ed array of ready flags for each thread */
        WMB();

        /* setup and create threads */
        for (i = 0; i < threads; i++) {
            ti[i].tid = i;
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
    free(quit_flag);

    cleanup_ctx();

    MPI_Comm_free(&comm);

    MPI_Finalize();

    return 0;
}
