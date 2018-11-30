/*
 * Copyright (c) 2017-2018 Mellanox Technologies Ltd. ALL RIGHTS RESERVED.
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#ifndef TIMELINE_H
#define TIMELINE_H

typedef enum {
    EVENT_POST,
    EVENT_WAIT
} mtcomb_event_t;

#define COLOR_1 "#00A2DE"
#define COLOR_2 "#101073"

typedef struct {
    mtcomb_event_t type;
    double start, end;
} timeline_event_t;

typedef struct {
    int thr_id;
    int num_events;
    timeline_event_t events[100];
} thread_timeline_t;


typedef struct {
    int proc_id;
    int thr_num;
    thread_timeline_t threads[8];
} proc_timeline_t;

typedef struct {
    int proc_num;
    proc_timeline_t procs[8];
} timeline_t;

typedef struct {
    int thr_id;
    int proc_id;
} resource_map_t;


#endif // TIMELINE_H
