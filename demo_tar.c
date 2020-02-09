#define STDARC_C
#include "stdarc.c"

#include <stdio.h>

int main(int argc, char **argv) {
    const char *fname = argc > 1 ? argv[1] : "demo.tar";
    printf("list contents of %s ...\n", fname);
    tar *t = tar_open(fname, "rb");
    if( t ) {
        for( unsigned i = 0 ; i < tar_count(t); ++i ) {
            printf("  %d) @%08x %11u %s ", i+1, tar_offset(t,i), tar_size(t,i), tar_name(t,i));
            void *data = tar_extract(t,i);
            printf("\r%c\n", data ? 'Y':'N'); // use data here: "%.*s\n",tar_size(t,i),(char*)data
            if(data) free(data); 
        }
        tar_close(t);
    }
}
