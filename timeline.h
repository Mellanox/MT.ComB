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

#include <stdio.h>

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
#ifdef STANDALONE
    timeline_event_t events[100];
#else
    timeline_event_t *events;
#endif
} thread_timeline_t;

typedef struct {
    int proc_id;
    int thr_num;
#ifdef STANDALONE
    thread_timeline_t threads[8];
#else
    thread_timeline_t *threads;
#endif
} proc_timeline_t;

typedef struct {
    int proc_num;
#ifdef STANDALONE
    proc_timeline_t procs[8];
#else
    proc_timeline_t *procs;
#endif
} timeline_t;

typedef struct {
    int thr_id;
    int proc_id;
} resource_map_t;

int write_timeline(timeline_t *tl, char *fname,
                   int start_p, int start_thr, int start_ev, int ev_count);
int serialize_timeline(timeline_t *tl, char *fname);
int deserialize_timeline(timeline_t *tl, char *fname);

#endif // TIMELINE_H
