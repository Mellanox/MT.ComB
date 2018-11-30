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
}

void special_usage(char *cmd) {
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
}

void *init_rreqs(int num_rreqs, void *data_buf, int msg_sz, int tag) {
}

void cleanup_sreqs(void *sreqs) {
}

void cleanup_rreqs(void *rreqs) {
}


void nonblocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *sreqs, int idx) {
}

void nonblocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag, void *rreqs, int idx) {
}

void blocking_send(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag){

}
void blocking_recv(void* data_buf, int msg_sz, struct thread_info *tinfo, int tag){
}

void wait_all_sreqs(void *sreqs, struct thread_info *tinfo, int num_sreqs){
}

void wait_all_rreqs(void *rreqs, struct thread_info *tinfo, int num_rreqs){
}

int main(int argc, char *argv[]) {
}
