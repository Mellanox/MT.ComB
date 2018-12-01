#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "timeline.h"

int procnum = -1, thrnum = -1, evnum = -1, evcnt = -1 /* max possible */;
char *input = NULL, *output = NULL;

void usage(char *cmd) {
    fprintf(stderr, "Options: %s\n", cmd);
    fprintf(stderr, "\t-i\tInput file name\n");
    fprintf(stderr, "\t-o\tOutput file name\n");
    fprintf(stderr, "\t-p <procnum>\tStart from the event belonging to process <procnum>\n");
    fprintf(stderr, "\t-t <thrnum>\tStart from the event belonging to thread <thrnum>\n");
    fprintf(stderr, "\t-e <evnum>\tStart from the event number <evnum>\n");
    fprintf(stderr, "\t-n <evcnt>\tDisplay <evcnt> events on monitored thread\n");
}

int check_unsigned(char *str) {
    return (strlen(str) == strspn(str,"0123456789") );
}

void process_args(int argc, char **argv) {
    int c, rank, rc = 0;

    while ((c = getopt(argc, argv, "hi:o:p:t:e:n:")) != -1) {
        switch (c) {
        case 'h':
                usage(argv[0]);
            break;
        case 'p':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            procnum = atoi(optarg);
            break;
        case 't':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            thrnum = atoi(optarg);
            break;
        case 'e':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            evnum = atoi(optarg);
            break;
        case 'n':
            if (!check_unsigned(optarg)) {
                goto error;
            }
            evcnt = atoi(optarg);
            break;
        case 'i':
            input = strdup(optarg);
            break;
        case 'o':
            output = strdup(optarg);
            break;
        default:
            c = -1;
            goto error;
        }
    }
    return;
error:
    if(c != -1) {
        fprintf(stderr, "Bad argument of '-%c' option\n", (char)c);
    }
    usage(argv[0]);
    exit(1);
}


int main(int argc, char **argv)
{
    timeline_t tl;

    process_args(argc, argv);
    if( input == NULL || output == NULL ){
        printf("Need input and output files\n");
        exit(1);
    }

    deserialize_timeline(&tl, input);
    write_timeline(&tl, output, procnum, thrnum, evnum, evcnt);
}
