/*
 * Copyright (c) 2017-2018 Mellanox Technologies Ltd. ALL RIGHTS RESERVED.
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

//#define STANDALONE
#include "timeline.h"

int ystep = 1;
float ydelta = 0.4;

//#define STANDALONE
#ifdef STANDALONE
timeline_t tl = {
    1,
    { /* procs */
      { /* proc #1 */
        0, /* proc id */
        2, /* thread num */
        { /* threads */
          { /* thread #1 */
              0, 4, /* 4 events */
              { /* events */
                  { EVENT_POST, 0, 0.1 },
                  { EVENT_WAIT, 0.11, 0.3 },
                  { EVENT_POST, 0.31, 0.5 },
                  { EVENT_WAIT, 0.51, 0.9 }
              }
          },
          { /* thread #2 */
              1, 4, /* 4 events */
              { /* events */
                  { EVENT_POST, 0, 0.11 },
                  { EVENT_WAIT, 0.12, 0.6 },
                  { EVENT_POST, 0.61, 0.8 },
                  { EVENT_WAIT, 0.81, 1.2 }
              }
          }
        }
      }
    }
};
#endif

void find_thread_xrange(thread_timeline_t *thr_tl, double *start, double *end)
{
    *start = thr_tl->events[0].start;
    *end = thr_tl->events[thr_tl->num_events-1].end;
}

void find_proc_xrange(proc_timeline_t *proc_tl, double *_start, double *_end)
{
    int t;
    double start = -1, end = -1;

    for(t = 0; t < proc_tl->thr_num; t++){
        double t_start, t_end;
        find_thread_xrange(&proc_tl->threads[t], &t_start, &t_end);
        if( 0 > start || (start < t_start) ) {
            start = t_start;
        }
        if( 0 > end || (end < t_end) ) {
            end = t_end;
        }
    }
    *_start = start;
    *_end = end;
}

int find_timeline_xrange(timeline_t *tl, double *_start, double *_end,
                          int ref_proc, int ref_thr, int ref_ev, int ev_cnt)
{
    int p;
    double start = -1, end = -1;

    if( (ref_proc > tl->proc_num) || (ref_thr > tl->procs[ref_proc].thr_num)
            || (ref_ev > tl->procs[ref_proc].threads[ref_thr].num_events)
            || ((ref_ev + ev_cnt) > tl->procs[ref_proc].threads[ref_thr].num_events) ){
        return EINVAL;
    }

    if( (0 > ref_proc) || (0 > ref_thr) || (0 > ref_ev) || (0 > ev_cnt)) {
        for(p = ref_proc; p < tl->proc_num; p++){
            double t_start, t_end;
            find_proc_xrange(&tl->procs[p], &t_start, &t_end);
            if( 0 > start || (start > t_start) ) {
                start = t_start;
            }
            if( 0 > end || (end < t_end) ) {
                end = t_end;
            }
        }
    } else {
        start = tl->procs[ref_proc].threads[ref_thr].events[ref_ev].start;
        end = start = tl->procs[ref_proc].threads[ref_thr].events[ref_ev + ev_cnt].end;
    }
    *_start = start;
    *_end = end;
}

int find_timeline_res(timeline_t *tl, resource_map_t **_map)
{
    int p, t, nres = 0, res_idx = 0;
    resource_map_t *map = NULL;

    for(p = 0; p < tl->proc_num; p++){
        nres += tl->procs[p].thr_num;
    }

    map = calloc(nres, sizeof(*map));

    for(p = 0; p < tl->proc_num; p++){
        for(t = 0; t < tl->procs[p].thr_num; t++) {
            map[res_idx].proc_id = p;
            map[res_idx++].thr_id = t;
        }
    }
    *_map = map;
    return nres;
}

void build_yticks(timeline_t *tl, resource_map_t *map, int nres, int step,
                 char *ytics, size_t len)
{
    int i = 0;
    char *ptr = ytics;
    int no_proc = 0;
    *(ptr++) = '(';
    *(ptr) = '\0';

    if( 1 == tl->proc_num ) {
        no_proc = 1;
    }

    for(i = 0; i < nres; i++){
        if( no_proc ) {
            ptr += snprintf(ptr, len - (ptr - ytics), "\"thr%d\" %d",
                            map[i].thr_id, step * (i + 1) );
        } else {
            ptr += snprintf(ptr, len - (ptr - ytics), "\"p%d.t%d\" %d",
                            map[i].proc_id, map[i].thr_id,
                            step * (i + 1) );
        }
        if( i < (nres - 1) ){
            ptr += snprintf(ptr, len - (ptr - ytics), ", ");
        }
    }

    snprintf(ptr, len - (ptr - ytics), " )");
}


void write_bars(FILE *fp, timeline_t *tl, double r_start, double r_end,
                resource_map_t *map, int nres, int ystep)
{
    int p, t;
    int objnum = 1;

    for(p=0; p < tl->proc_num; p++){
        for(t=0; t < tl->procs[p].thr_num; t++){
            int i = -1, e;
            for(i = 0; i < nres; i++){
                if( map[i].proc_id == p && map[i].thr_id == t ){
                    break;
                }
            }
            for(e = 0; e < tl->procs[p].threads[t].num_events; e++){
                timeline_event_t *ev = &tl->procs[p].threads[t].events[e];
                /* Filter by the time range */
                if( ev->start < r_start || ev->end < r_end ){
                    continue;
                }
                fprintf(fp, "set object %d rectangle from ", objnum++);
                fprintf(fp, "%.3lf, %f",
                        1E6*(ev->start - r_start),
                        (i + 1) * ystep - ydelta);
                fprintf(fp," to ");
                fprintf(fp, "%.3lf, %f",
                        1E6*(ev->end - r_start),
                        (i + 1) * ystep + ydelta);
                fprintf(fp, " fillcolor rgb ");
                switch(ev->type) {
                case EVENT_POST:
                    fprintf(fp, "\"%s\"", COLOR_1);
                    break;
                case EVENT_WAIT:
                    fprintf(fp, "\"%s\"", COLOR_2);
                    break;
                }
                fprintf(fp, " fillstyle solid 0.8\n");
            }
        }
    }
}

char *prefix =
        "set terminal postscript eps color solid\n"
        "set output \"test.eps\"\n"
        "set autoscale x\n"
        "set xlabel \"time\"\n"
        "set ylabel \"\"\n"
        "set key outside width +2\n"
        "set grid xtics\n"
        "set palette model RGB\n"
        "unset colorbox\n";

int _write_timeline(timeline_t *tl, FILE *fp,
                    int start_p, int start_thr, int start_ev, int ev_count)
{
    double start, end;
    int nres;
    resource_map_t *map;

    /* Printf fixed prefix */
    fprintf(fp, "%s", prefix);

    if( (0 > start_p) || (0 > start_thr) || (0 > start_ev)){
        start_p = 0;
        start_thr = 0;
        start_ev = 0;
    }

    /* Print ranges and y lables */
    find_timeline_xrange(tl, &start, &end,
                         start_p, start_thr, start_ev, ev_count);
    fprintf(fp, "set xrange [0.000:%lf]\n", 1E6*(end - start));
    char ytics[1024] = "";
    nres = find_timeline_res(tl, &map);
    build_yticks(tl, map, nres, ystep, ytics, 1024);
    fprintf(fp, "set yrange [0.4:%lf]\n", (double)nres + 0.6);
    fprintf(fp, "set ytics %s\n", ytics);
    fprintf(fp, "set title \"MT.ComB: %d-thread timeline\"\n", nres);

    write_bars(fp, tl, start, end, map, nres, ystep);

    fprintf(fp, "plot \\\n"
            "-1 title \"Post\" with lines linecolor rgb \"" COLOR_1 "\" linewidth 6, \\\n"
            "-1 title \"Wait\" with lines linecolor rgb \"" COLOR_2 "\" linewidth 6\n");

}

int write_timeline(timeline_t *tl, char *fname,
                   int ref_proc, int ref_thr, int ref_ev, int ev_cnt)
{
    FILE *fp = fopen(fname, "w");
    _write_timeline(tl, fp, ref_proc, ref_thr, ref_ev, ev_cnt);
}


int serialize_item(void *ptr, int elem_size, int elem_count, FILE *fp)
{
    size_t written = 0;
    written = fwrite(ptr, elem_size, elem_count, fp);
    if( written != (elem_count * elem_size) ) {
        return EIO;
    }
    return 0;
}

int deserialize_item(void *ptr, int elem_size, int elem_count, FILE *fp)
{
    size_t read = 0;
    read = fread(ptr, elem_size, elem_count, fp);
    if( read != (elem_count * elem_size) ) {
        return EIO;
    }
    return 0;
}

#define SERIALIZE_ITEM(ptr, esize, ecnt, fp) {    \
    int rc;                                       \
    rc = serialize_item(ptr, esize, ecnt, fp);    \
    if( rc ){                                     \
        return rc;                                \
    }                                             \
}

#define DESERIALIZE_ITEM(ptr, esize, ecnt, fp) {    \
    int rc;                                         \
    rc = deserialize_item(ptr, esize, ecnt, fp);    \
    if( rc ){                                       \
        return rc;                                  \
    }                                               \
}


int serialize_thread(thread_timeline_t *ttl, FILE *fp)
{
    SERIALIZE_ITEM(&ttl->thr_id, sizeof(ttl->thr_id), 1, fp);
    SERIALIZE_ITEM(&ttl->num_events, sizeof(ttl->num_events), 1, fp);
    SERIALIZE_ITEM(&ttl->events, sizeof(*ttl->events), ttl->num_events, fp);
    return 0;
}

int deserialize_thread(thread_timeline_t *ttl, FILE *fp)
{
    DESERIALIZE_ITEM(&ttl->thr_id, sizeof(ttl->thr_id), 1, fp);
    DESERIALIZE_ITEM(&ttl->num_events, sizeof(ttl->num_events), 1, fp);
    if(!(ttl->events = calloc(ttl->num_events,sizeof(*ttl->events)))) {
        return -ENOMEM;
    }
    DESERIALIZE_ITEM(&ttl->events, sizeof(*ttl->events), ttl->num_events, fp);
    return 0;
}

int serialize_process(proc_timeline_t *ptl, FILE *fp)
{
    int i, rc;
    SERIALIZE_ITEM(&ptl->proc_id, sizeof(ptl->proc_id), 1, fp);
    SERIALIZE_ITEM(&ptl->thr_num, sizeof(ptl->thr_num), 1, fp);
    for(i=0; i<ptl->thr_num; i++){
        if(rc = serialize_thread(&ptl->threads[i], fp) ) {
            return rc;
        }
    }
    return 0;
}

int deserialize_process(proc_timeline_t *ptl, FILE *fp)
{
    int i, rc;
    DESERIALIZE_ITEM(&ptl->proc_id, sizeof(ptl->proc_id), 1, fp);
    DESERIALIZE_ITEM(&ptl->thr_num, sizeof(ptl->thr_num), 1, fp);
    ptl->threads = calloc(ptl->thr_num, sizeof(*ptl->threads));
    if(!(ptl->threads = calloc(ptl->thr_num,sizeof(*ptl->threads)))) {
        return -ENOMEM;
    }
    for(i=0; i<ptl->thr_num; i++){
        if(rc = deserialize_thread(&ptl->threads[i], fp) ) {
            return rc;
        }
    }
    return 0;
}


int serialize_timeline(timeline_t *tl, char *fname)
{
    FILE *fp = fopen(fname, "w");
    int i, rc;
    if( !fp ) {
        return EIO;
    }
    SERIALIZE_ITEM(&tl->proc_num, sizeof(tl->proc_num), 1, fp);
    for(i = 0; i < tl->proc_num; i++) {
        if( rc = serialize_process(&tl->procs[i], fp) ) {
            return rc;
        }
    }
    fclose(fp);
}

int deserialize_timeline(timeline_t *tl, char *fname)
{
    FILE *fp = fopen(fname, "r");
    int i, rc;
    if( !fp ) {
        return EIO;
    }
    SERIALIZE_ITEM(&tl->proc_num, sizeof(tl->proc_num), 1, fp);
    tl->procs = calloc(tl->proc_num, sizeof(*tl->procs));
    for(i = 0; i < tl->proc_num; i++) {
        if( rc = deserialize_process(&tl->procs[i], fp) ) {
            return rc;
        }
    }
    fclose(fp);
}


#ifdef STANDALONE
int main()
{
    write_timeline(&tl, "output.gpl", 0, 0, 0, -1);
}
#endif
