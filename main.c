#include "generic.h"

void set_default_args() {
    worker = worker_nb;
    threads = 1;
    win_size = 256;
    iterations = 50;
    warmup = 10;
    msg_size = 0;
}

int main(int argc, char *argv[]) {
    int rank, size, i;
    pthread_t *id;
    MPI_Comm comm;
    int ret, mt_level_act;

    set_default_args();

    pre_scan_args(argc, argv);

    if (want_thr_support && is_mpi_based_test()) {
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mt_level_act);
        if( mt_level_act != MPI_THREAD_MULTIPLE ){
            if (rank == 0) {
                fprintf(stderr, "ERROR: Unable to initialize MPI thread support!\n");
                MPI_Finalize();
                exit(0);
            }
        }
    } else {
        MPI_Init(&argc, &argv);
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

    if (ret = mem_map()) {
        MPI_Finalize();
        exit(1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (threads == 1) {
        setup_thread_info_single(ti);
        sync_start_all = sync_cur_step;
        worker((void*)ti);
    } else {
        /* Create the zero'ed array of ready flags for each thread */
        WMB();
        /* setup and create threads */
        for (i = 0; i < threads; i++) {
            setup_thread_info_multi(ti, i);
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

    cleanup_thread_info(ti, threads);

    free(id);
    free(results);
    free(global_buf);
    free(ti);
    cleanup_ctx();
    MPI_Comm_free(&comm);
    MPI_Finalize();
    return 0;
}
