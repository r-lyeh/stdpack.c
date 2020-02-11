#include <stdio.h>
#include <string.h>

#define STDARC_C
#include "stdarc.c"

int main(int argc, char **argv) {
    // append file to a zip. will create if does not exist
    puts("appending file to demo.pak ...");
    pak *z = pak_open("demo.pak", "a+b");
    if( z ) {
        const char *fname = __FILE__;
        for( FILE *fp = fopen(fname, "rb"); fp; fclose(fp), fp = 0) {
            pak_append_file(z, fname, fp);
        }
        pak_close(z);
    }

    // test contents of file
    char *fname = argc > 1 ? argv[1] : "demo.pak";
    printf("testing files in %s ...\n", fname);
    z = pak_open(fname, "rb");
    if( z ) {
        for( unsigned i = 0 ; i < pak_count(z); ++i ) {
            printf("  %d) ", i+1);
            printf("@%08x ", pak_offset(z,i));
            printf("%11u ", pak_size(z,i));
            printf("%s ", pak_name(z,i));
            void *data = pak_extract(z,i);
            // use data [...]
            printf("\r%c\n", data ? 'Y':'N'); // %.*s\n", pak_size(z,i), (char*)data);
            if(data) free(data); 
        }
        pak_close(z);
    }

    puts("ok");
}
