/*
 * Copyright (c) 2017 Mellanox Technologies Ltd. ALL RIGHTS RESERVED.
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "generic.h"

int threads;
int win_size;
int iterations;
int warmup;
int dup_comm;
int msg_size;
int binding = 0;
int want_thr_support = MODE_MT_BASIC;
int intra_node = 0;
int verify_mode = 0;
int send_alignment = 1;
int recv_alignment = 1;
int buf_unit_size, remote_buf_unit_size;

int my_host_idx = -1, my_rank_idx = -1, my_leader = -1;
int my_partner, i_am_sender = 0;

/* binding-related */
cpu_set_t my_cpuset;
long configured_cpus;

/* Message rate for each thread */
double *results;

/*data buf on sendera and receiver */
char *global_buf = NULL;

/* thread synchronization */
int volatile sync_cur_step = 1;
int volatile *sync_thread_ready = NULL;
int volatile sync_start_all = 0;

char *binding_err_msgs[] = {
    "all is OK",
    "some nodes have more procs than CPUs. Oversubscription is not supported by now",
    "some process have overlapping but not equal bindings. Not supported by now",
    "some process was unable to read it's affinity",
};

void sync_worker(struct thread_info *tinfo) {
    int tid = tinfo->tid;
    if (threads == 1) {
        wireup_progress(tinfo);
        MPI_Barrier(MPI_COMM_WORLD);
    } else {
        int cur_step = sync_cur_step;

        WMB();

        /* now we are ready to send */
        sync_thread_ready[tid] = sync_cur_step;

        WMB();

        /* wait for master's signal */
        while (sync_start_all != cur_step) {
            usleep(1);
            wireup_progress(tinfo);
        }
    }
}

void sync_master() {
    int i;
    int cur_step = sync_cur_step;

    /* wait for all threads to start */
    WMB();

    for (i = 0; i < threads; i++) {
        while (sync_thread_ready[i] != cur_step) {
            usleep(10);
        }
    }

    /* Synchronize all processes */
    WMB();
    MPI_Barrier(MPI_COMM_WORLD);
    sync_cur_step++;

    /* signal threads to start */
    WMB();
    sync_start_all = cur_step;
}

void pre_scan_args(int argc, char **argv) {
    int i;
    /* this is a hack - you cannot access argc/argv before a call to MPI_Init[_thread].
     * But it seems to work with most modern MPIs: Open MPI and Intel MPI was checked
     */
    for (i = 0; i < argc; i++) {
        if(0 == strcmp(argv[i], "-Dthrds")) {
            want_thr_support = MODE_MT_SINGLE;
        }
    }
}

void usage(char *cmd) {
    set_default_args();

    fprintf(stderr, "Options: %s\n", cmd);
    fprintf(stderr, "\t-h\tDisplay this help\n");
    fprintf(stderr, "Test description:\n");
    fprintf(stderr, "\t-s\tMessage size (default: %d)\n", msg_size);
    fprintf(stderr, "\t-n\tNumber of measured iterations (default: %d)\n", iterations);
    fprintf(stderr, "\t-w\tNumber of warmup iterations (default: %d)\n", warmup);
    fprintf(stderr, "\t-W\tWindow size - number of send/recvs between sync (default: %d)\n", win_size);
    fprintf(stderr, "\t-t\tNumber of threads (default: %d)\n", threads);
    fprintf(stderr, "\t-l\tAlignment on sender buffer (default: %d\n", send_alignment);
    fprintf(stderr, "\t-r\tAlignment on receiver buffer (default: %d\n", recv_alignment);
    fprintf(stderr, "Test options:\n");
    fprintf(stderr, "\t-B\tBlocking send (default: %s)\n", (worker == worker_b) ? "blocking" : "non-blocking");
    fprintf(stderr, "\t-S\tSMP mode - intra-node performance - pairwise exchanges (default: %s)\n", 
            (intra_node) ? "enabled" : "disabled" );
    fprintf(stderr, "\t-b\tEnable fine-grained binding of the threads inside the set provided by MPI\n");
    fprintf(stderr, "\t\tbenchmark is able to discover node-local ranks and exchange existing binding\n");
    fprintf(stderr, "\t\t(default: %s)\n", (binding) ? "enabled" : "disabled");
    fprintf(stderr, "\t-v\tEnable verification mode (default %d)\n", verify_mode);
    special_usage(cmd);
}


int check_unsigned(char *str) {
    return (strlen(str) == strspn(str,"0123456789") );
}


void process_args(int argc, char **argv) {
    int c, rank, rc = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    while ((c = getopt(argc, argv, "hD:t:BbW:n:w:ds:Svl:r:E:")) != -1) {
        switch (c) {
        case 'h':
            if (0 == rank) {
                usage(argv[0]);
            }
            MPI_Finalize();
            exit(0);
            break;
        case 'D':
            /* from the getopt perspective thrds is an optarg. Check that 
             * user haven't specified anything else */
            if (strcmp(optarg, "thrds")) {
                goto error;
            }
            break;
        case 't':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            threads = atoi(optarg);
            if (threads == 0) {
                goto error;
            }
            break;
        case 'B':
            worker = worker_b;
            break;
        case 'b':
            binding = 1;
            break;
        case 'W':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            win_size = atoi(optarg);
            if (win_size == 0) {
                goto error;
            }
            break;
        case 'n':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            iterations = atoi(optarg);
            if (iterations == 0) {
                goto error;
            }
            break;
        case 'w':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            warmup = atoi(optarg);
            if (warmup < 0) {
                goto error;
            }
            break;
        case 's':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            msg_size = atoi(optarg);
            if (msg_size == 0) {
                goto error;
            }
            break;
        case 'd':
            dup_comm = 1;
            break;
        case 'S':
            intra_node = 1;
            break;
        case 'v':
            verify_mode = 1;
            break;
        case 'l':
            send_alignment = atoi(optarg);
            if (send_alignment < 1) {
                goto error;
            }
            break;
        case 'r':
            recv_alignment = atoi(optarg);
            if (recv_alignment < 1) {
                goto error;
            }
            break;
        case 'E':
            if (special_process_args(optarg)) {
                goto error;
            }
            break;
        default:
            c = -1;
            goto error;
        }
    }
    return;
error:
    if(0 == rank) {
        if(c != -1) {
            fprintf(stderr, "Bad argument of '-%c' option\n", (char)c);
        }
        usage(argv[0]);
    }
    MPI_Finalize();
    exit(1);
}

void setup_binding(MPI_Comm comm) {
    int grank;
    cpu_set_t set, *all_sets;
    int *ranks, ranks_cnt = 0, my_idx = -1;
    int rank, size, ncpus;
    int error_loc = 0, error;
    int cpu_no, cpus_used = 0;
    int i;

    if (!binding) {
        return;
    }

    /* How many cpus we have on the node */
    configured_cpus = sysconf(_SC_NPROCESSORS_CONF);

    /* MPI identifications */
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    ranks = calloc(size, sizeof(*ranks));

    if (sizeof(set) * 8 < configured_cpus) {
        /* Shouldn't happen soon, currentlu # of cpus = 1024 */
        if (0 == grank) {
            fprintf(stderr, "The size of CPUSET is larger that we can currently handle\n");
        }
        MPI_Finalize();
        exit(1);
    }

    if (sched_getaffinity(0, sizeof(set), &set)) {
        error_loc = BIND_GETAFF;
    }

    ncpus = CPU_COUNT(&set);

    /* FIXME: not portable, do we care about heterogenious systems? */
    all_sets = calloc(size, sizeof(set));
    MPI_Allgather(&set, sizeof(set), MPI_BYTE, all_sets, sizeof(set), MPI_BYTE, comm);

    if (error_loc) {
        goto finish;
    }

    /* find procs that are bound exactly like we are.
     * Note that we don't expect overlapping set's:
     *  - either they are the same;
     *  - or disjoint.
     */
    for (i = 0; i < size; i++) {
        cpu_set_t *rset = all_sets + i, cmp_set;
        CPU_AND(&cmp_set, &set, rset);
        if (CPU_COUNT(&cmp_set) == ncpus) {
            /* this rank is binded as we are */
            ranks[ranks_cnt++] = i;
        } else if (!CPU_COUNT(&cmp_set)) {
            /* other binding */
            continue;
        } else {
            /* not expected. Error exit! */
            error_loc = BIND_OVERLAP;
            goto finish;
        }
    }

    if (ncpus < ranks_cnt * threads) {
        error_loc = BIND_OVERSUBS;
    }

    for (i = 0; i < ranks_cnt; i++) {
        if (ranks[i] == rank) {
            my_idx = i;
            break;
        }
    }

    /* sanity check */
    assert(my_idx >= 0);

    cpu_no = 0;
    CPU_ZERO(&my_cpuset);
    for (i = 0; i < configured_cpus && cpus_used < threads; i++) {
        if (CPU_ISSET(i, &set)) {
            if (cpu_no >= (my_idx * threads)) {
                CPU_SET(i, &my_cpuset);
                cpus_used++;
            }
            cpu_no++;
        }
    }

finish:
    MPI_Allreduce(&error_loc, &error, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (error) {
        if (0 == grank) {
            fprintf(stderr, "ERROR (binding): %s\n", binding_err_msgs[error]);
        }
    }
}

void bind_worker(int tid) {
    int i, cpu_no = 0;

    if (!binding) {
        return;
    }

    for (i = 0; i < configured_cpus; i++) {
        if (CPU_ISSET(i, &my_cpuset)) {
            if (cpu_no == tid) {
                /* this core is mine! */
                cpu_set_t set;
                CPU_ZERO(&set);
                CPU_SET(i, &set);
                if (pthread_setaffinity_np(pthread_self(), sizeof(set), &set)) {
                    MPI_Abort(MPI_COMM_WORLD, 0);
                }
                break;
            }
            cpu_no++;
        }
    }

#ifdef DEBUG
    if (0 == my_host_idx) {
        int rank;
        cpu_set_t set;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (pthread_getaffinity_np(pthread_self(), sizeof(set), &set)) {
            MPI_Abort(MPI_COMM_WORLD, 0);
        }
        for (i = 0; i < configured_cpus; i++) {
            if (CPU_ISSET(i, &set)) {
                fprintf(stderr, "%d:%d: %d\n", rank, tid, i);
            }
        }
    }
#endif
}

/* Split ranks to the pairs communicating with each-other.
 * Following restrictions allpied (may be relaxed if need):
 * - there can only be 2 hosts, error otherwise
 * - number of ranks on each host must be equal
 * otherwise - err exit.
 */

MPI_Comm split_to_pairs() {
    int rank, size, len, max_len;
    char hname[MAX_HOSTHAME], *all_hosts, hnames[MAX_HOSTS][MAX_HOSTHAME];
    int i, j, hostnum = 0, *hranks[MAX_HOSTS], hranks_cnt[MAX_HOSTS] = {0, 0};
    MPI_Group base_grp, my_grp;
    MPI_Comm comm, bind_comm;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size % 2) {
        if (rank == 0) {
            fprintf(stderr,"Expect even number of ranks!\n");
        }
        MPI_Finalize();
        exit(1);
    }

    gethostname(hname, 1024);
    len = strlen(hname) + 1; /* account '\0' */
    MPI_Allreduce(&len, &max_len, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    all_hosts = calloc(sizeof(char), max_len * size);
    MPI_Allgather(hname, max_len, MPI_CHAR, all_hosts, max_len, MPI_CHAR, MPI_COMM_WORLD);

    for (i = 0; i < MAX_HOSTS; i++) {
        hranks[i] = calloc(sizeof(int), size);
    }

    for (i = 0; i < size; i++) {
        char *base = all_hosts + i * max_len;
        int hidx = -1;
        for (j = 0; j < hostnum; j++) {
            if (!strcmp(hnames[j], base)) {
                /* host already known, add rank */
                hidx = j;
                goto save_rank;
            }
        }
        /* new host */
        if (hostnum == 2) {
            /* too many hosts */
            if (rank == 0) {
                fprintf(stderr,"Please, launch your application on 2 hosts!\n");
            }
            MPI_Finalize();
            exit(1);
        }
        strcpy(hnames[hostnum], base);
        hidx = hostnum;
        hostnum++;
save_rank:
        hranks[hidx][hranks_cnt[hidx]] = i;
        if (i == rank) {
            my_host_idx = hidx;
            my_rank_idx = hranks_cnt[hidx];
        }
        hranks_cnt[hidx]++;
    }

    /* sanity check, this is already ensured by  */
    if (!intra_node)
        assert(hostnum == 2);
    else
        assert(hostnum == 1);

#ifdef DEBUG
    if (0 == rank) {
        /* output the rank-node mapping */
        for (i = 0; i < hostnum; i++) {
            int j;
            printf("%s: ", hnames[i]);
            for (j = 0; j < hranks_cnt[i]; j++) {
                printf("%d ", hranks[i][j]);
            }
            printf("\n");
        }
    }
#endif

    if (!intra_node) {
        /* sanity check */
        if (hranks_cnt[0] != (size / 2)) {
            if (rank == 0) {
                fprintf(stderr,"Ranks are non-evenly distributed on the nodes!\n");
            }
            MPI_Finalize();
            exit(1);
        }

        /* return my partner */
        i_am_sender = 0;
        if (my_host_idx == 0) {
            i_am_sender = 1;
        }
        my_partner = hranks[(my_host_idx+1)%MAX_HOSTS][my_rank_idx];
        my_leader = hranks[my_host_idx][0];

        /* create the communicator for all senders */
        MPI_Comm_group(MPI_COMM_WORLD, &base_grp);
        MPI_Group_incl(base_grp, size/2, hranks[my_host_idx], &my_grp);
        MPI_Comm_create(MPI_COMM_WORLD, my_grp, &comm);

        /* FIXME: do we need to free it here? Do we need to free base_grp? */
        MPI_Group_free(&my_grp);
        bind_comm = comm;
    } else {
        /* sanity check */
        if (hranks_cnt[0] != size) {
            if (rank == 0) {
                fprintf(stderr,"Ranks are non-evenly distributed on the nodes!\n");
            }
            MPI_Finalize();
            exit(1);
        }

        /* return my partner */
        i_am_sender = 0;
        if (my_rank_idx % 2 == 0) {
            i_am_sender = 1;
        }
        my_partner = hranks[0][my_rank_idx+(1-2*(my_rank_idx%2))];
        my_leader = my_rank_idx % 2;

        /* create the communicator for all senders */
        MPI_Comm_split(MPI_COMM_WORLD, my_rank_idx%2, my_rank_idx, &comm);
        bind_comm = MPI_COMM_WORLD;
    }

    /* discover and exchange binding info */
    setup_binding(bind_comm);

    /* release the resources */
    free(all_hosts);
    for (i = 0; i < hostnum; i++) {
        free(hranks[i]);
    }

    return comm;
}


void *worker_nb(void *info) {
    struct thread_info *tinfo = (struct thread_info*)info;
    int tag, i, j;
    double stime, etime;
    char *databuf = NULL;

    databuf = &(global_buf[tinfo->tid * buf_unit_size]);
    tag = tinfo->tid;
    bind_worker(tag);

    /* start the benchmark */
    if (i_am_sender) {
        void *sreqs = init_sreqs(win_size, (void *)databuf, msg_size, tag);

        for (i = 0; i < (iterations + warmup); i++) {
            if (i == warmup) {
                /* ensure that all threads start "almost" together */
                sync_worker(tinfo);
                stime = MPI_Wtime();
            }

            for (j = 0; j < win_size; j++) {
                nonblocking_send(databuf, msg_size, tinfo, tag, sreqs, j);
            }

            wait_all_sreqs(sreqs, tinfo, win_size);
        }

        cleanup_sreqs(sreqs);
    } else {
        void *rreqs = init_rreqs(win_size, (void *)databuf, msg_size, tag);

        for (i = 0; i < (iterations + warmup); i++) {
            if (i == warmup) {
                /* ensure that all threads start "almost" together */
                sync_worker(tinfo);
                stime = MPI_Wtime();
            }

            for (j = 0; j < win_size; j++) {
                nonblocking_recv(databuf, msg_size, tinfo, tag, rreqs, j);
            }

            wait_all_rreqs(rreqs, tinfo, win_size);
        }

        cleanup_rreqs(rreqs);
    }

    etime = MPI_Wtime();

    results[tag] = 1 / ((etime - stime) / (iterations * win_size));

    return 0;
}


/*
 * Derived from the threaded test written by R. Thakur and W. Gropp:
 * URL: http://www.mcs.anl.gov/~thakur/thread-tests/
 */
void *worker_b(void *info) {
    struct thread_info *tinfo = (struct thread_info*)info;
    int tag, i, j;
    double stime, etime;
    char *databuf = &(global_buf[tinfo->tid * buf_unit_size]);

    tag = tinfo->tid;
    bind_worker(tag);

    if (i_am_sender) {
        for (i = 0; i < iterations + warmup; i++) {
            if (i == warmup) {
                /* ensure that all threads start "almost" together */
                sync_worker(tinfo);
                stime = MPI_Wtime();
            }

            for (j = 0; j < win_size; j++) {
                blocking_send(databuf, msg_size, tinfo, tag);
            }
        }
    } else {
        for (i = 0; i < iterations + warmup; i++) {
            if (i == warmup) {
                /* ensure that all threads start "almost" together */
                sync_worker(tinfo);
                stime = MPI_Wtime();
            }

            for (j = 0; j < win_size; j++) {
                blocking_recv(databuf, msg_size, tinfo, tag);
            }
        }
    }

    etime = MPI_Wtime();
    results[tag] = 1 / ((etime - stime) / (iterations * win_size));

    return 0;
}

void verify_buf() {
    int i, k;
    for (i = 0; i < threads; i++) {
        for (k = 0; k < msg_size; k++) {
            if(global_buf[i*buf_unit_size+k] != i) {
                if (!i_am_sender) {
                    printf("ERROR: Recv buffer was corrupted at byte %d: "
                           "expect: 0x%hhx, found 0x%hhx\n",
                           i*buf_unit_size+k, i, global_buf[i*buf_unit_size+k]);
                } else {
                    printf("ERROR: Send buffer was corrupted at byte %d: "
                           "expect: 0x%hhx, found 0x%hhx\n",
                           i*buf_unit_size+k, i, global_buf[i*buf_unit_size+k]);
                }
            }
        }
    }
}

void allocate_global_buf() {
    int align, i, j;
    if (i_am_sender) {
        align = send_alignment;
        remote_buf_unit_size = (msg_size / recv_alignment + !!(msg_size % recv_alignment)) * recv_alignment;
    } else {
        align = recv_alignment;
    }

    buf_unit_size = (msg_size / align + !!(msg_size % align)) * align;
    global_buf = aligned_alloc(buf_unit_size, threads * buf_unit_size);

    if (verify_mode) {
        if (i_am_sender) {
            for (i = 0; i < threads; i++) {
                for (j = 0; j < msg_size; j++) {
                    global_buf[i*buf_unit_size + j] = i;
                }
            }
        } else {
            memset(global_buf, 0, threads * buf_unit_size);
        }
    }
}

void print_results(MPI_Comm comm) {
#ifdef DEBUG
    char tmp[1024] = "";
    sprintf(tmp, "%d: ", rank);
    for (i = 0; i < threads; i++) {
        sprintf(tmp, "%s%lf ", tmp, results[i]);
    }
    printf("%s\n", tmp);
#endif

    if (i_am_sender) {
        /* FIXME: for now only count on sender side, extend if makes sense */
        double results_rank = 0, results_node = 0;
        int i;

        for (i = 0; i < threads; i++) {
            results_rank += results[i];
        }

        MPI_Reduce(&results_rank, &results_node, 1, MPI_DOUBLE, MPI_SUM, my_leader, comm);

        if (my_rank_idx == 0) { /* output: msg throughput per rank */
            printf(">\t%d\t%lf\n", msg_size, results_node);
        }
    }
}
