#include <stdio.h>
#include "getopt.h"


/*LINTLIBRARY*/

#ifndef NULL
#define NULL    0
#endif

extern short strcmp();
extern char *strchr();

INT_4   opterr = 1;
INT_4   optind = 1;
INT_4   optopt;
char    *optarg;

INT_4 getopt(INT_4 argc, CHAR **argv, CHAR *opts)
{
        static INT_4 sp = 1;
        register INT_4 c;
        register CHAR *cp;

        if(sp == 1)
                if(optind >= argc ||
                   argv[optind][0] != '-' || argv[optind][1] == '\0')
                        return(EOF);
                else if(strcmp(argv[optind], "--") == NULL) {
                        optind++;
                        return(EOF);
                }
        optopt = c = argv[optind][sp];
        if(c == ':' || (cp=strchr(opts, c)) == NULL) {
            fprintf(stderr,"illegal command-line option: %c\n",c);
                if(argv[optind][++sp] == '\0') {
                        optind++;
                        sp = 1;
                }
                return('?');
        }
        if(*++cp == ':') {
                if(argv[optind][sp+1] != '\0')
                        optarg = &argv[optind++][sp+1];
                else if(++optind >= argc) {
                    fprintf(stderr,"command-line option %c requires an argument\n",c);
                        sp = 1;
                        return('?');
                } else
                        optarg = argv[optind++];
                sp = 1;
        } else {
                if(argv[optind][++sp] == '\0') {
                        sp = 1;
                        optind++;
                }
                optarg = NULL;
        }
        return(c);
}
