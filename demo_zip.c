#include <stdio.h>
#include <string.h>

#define STDARC_C
#include "stdarc.c"

int main(int argc, char **argv) {
    // append file to a zip. will create if does not exist
    puts("appending file to demo.zip ...");
    zip *z = zip_open("demo.zip", "a+b");
    if( z ) {
        // compress with DEFLATE|6. Other compressors are also supported (try LZMA|5, ULZ|9, LZ4X|3, etc.)
        for( FILE *myfile = fopen(__FILE__, "rb"); myfile; fclose(myfile), myfile = 0) {
            zip_append_file(z, __FILE__, myfile, 6);
        }
        zip_close(z);
    }

    // test contents of file
    puts("testing files in demo.zip ...");
    z = zip_open("demo.zip", "rb");
    if( z ) {
        for( unsigned i = 0 ; i < zip_count(z); ++i ) {
            printf("  %d) ", i+1);
            printf("[%08X] ", zip_hash(z,i));
            printf("$%02X ", zip_codec(z,i));
            printf("%s ", zip_modt(z,i));
            printf("%5s ", zip_file(z,i) ? "" : "<dir>");
            printf("%11u ", zip_size(z,i));
            printf("%s ", zip_name(z,i));
            //printf("@%x ", zip_offset(z,i));
            void *data = zip_extract(z,i);
            // use data [...]
            printf("\r%c\n", data ? 'Y':'N'); // %.*s\n", zip_size(z,i), (char*)data);
            if(data) free(data); 
        }
        zip_close(z);
    }

    puts("ok");
}
