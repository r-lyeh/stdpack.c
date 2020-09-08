#include <stdio.h>
#include <string.h>

#define STDPACK_C
#include "../stdpack.c"

int main(int argc, char **argv) {
    const char *src = "Hello world! Hello world! Hello world! Hello world!";
    unsigned srclen = strlen(src)+1;

    char out[128];
    unsigned outlen = mem_encode(src, srclen, out, 128, LZ4X|6); // use LZ4X level 6
    printf("%s %u->%u\n", outlen ? "ok" : "fail", srclen, outlen);

    char redo[128];
    unsigned unpacked = mem_decode(out, outlen, redo, 128);
    printf("%u->%u %s\n", outlen, unpacked, redo);
}
