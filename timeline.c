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

void find_timeline_xrange(timeline_t *tl, double *_start, double *_end)
{
    int p;
    double start = -1, end = -1;

    for(p = 0; p < tl->proc_num; p++){
        double t_start, t_end;
        find_proc_xrange(&tl->procs[p], &t_start, &t_end);
        if( 0 > start || (start > t_start) ) {
            start = t_start;
        }
        if( 0 > end || (end < t_end) ) {
            end = t_end;
        }
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


void write_bars(FILE *fp, timeline_t *tl, double range_start,
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
                fprintf(fp, "set object %d rectangle from ", objnum++);
                fprintf(fp, "%.3lf, %f",
                        1E6*(ev->start - range_start),
                        (i + 1) * ystep - ydelta);
                fprintf(fp," to ");
                fprintf(fp, "%.3lf, %f",
                        1E6*(ev->end - range_start),
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

int _write_timeline(timeline_t *tl, FILE *fp)
{
    double start, end;
    int nres;
    resource_map_t *map;

    /* Printf fixed prefix */
    fprintf(fp, "%s", prefix);

    /* Print ranges and y lables */
    find_timeline_xrange(tl, &start, &end);
    fprintf(fp, "set xrange [0.000:%lf]\n", 1E6*(end - start));
    char ytics[1024] = "";
    nres = find_timeline_res(tl, &map);
    build_yticks(tl, map, nres, ystep, ytics, 1024);
    fprintf(fp, "set yrange [0.4:%lf]\n", (double)nres + 0.6);
    fprintf(fp, "set ytics %s\n", ytics);
    fprintf(fp, "set title \"MT.ComB: %d-thread timeline\"\n", nres);

    write_bars(fp, tl, start, map, nres, ystep);

    fprintf(fp, "plot \\\n"
            "-1 title \"Post\" with lines linecolor rgb \"" COLOR_1 "\" linewidth 6, \\\n"
            "-1 title \"Wait\" with lines linecolor rgb \"" COLOR_2 "\" linewidth 6\n");

}

int write_timeline(timeline_t *tl, char *fname)
{
    FILE *fp = fopen(fname, "w");
    _write_timeline(tl, fp);
}

#ifdef STANDALONE
int main()
{
    write_timeline(&tl, "output.gpl");
}
#endif
