#define STDARC_C
#include "stdarc.c"

#include <stdio.h>

int main(int argc, char **argv) {
    const char *fname = argc > 1 ? argv[1] : "./";
    printf("list contents of %s ...\n", fname);
    dir *t = dir_open(fname, "rb");
    if( t ) {
        for( unsigned i = 0 ; i < dir_count(t); ++i ) {
            printf("Y %3d) %11u %s\t%s\n", i+1, dir_size(t,i), dir_name(t,i), dir_file(t,i) ? "" : "<dir>");
        }
        dir_close(t);
    }
}
