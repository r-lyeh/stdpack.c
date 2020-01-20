int main() {
    include("begin.c");
        include("balz.c");
        include("bcm_bwt.c");
        include("bcm.c");
        include("crush.c");
        include("deflate.c");
        include("lz4x.c");
        include("lzma.c");
        include("lzp1.c");
        include("lzrw3a.c");
        include("lzss.c");
        include("ppp.c");
        include("raw.c");
        include("ulz.c");
    include("end.c");
}

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

int include( const char *file ) {
    struct stat st = {0};
    if( stat(file, &st) < 0 ) {
        fprintf(stderr, "Warning: Cannot open '%s'\n", file);
        return -1;
    }
    char cmd[256] = {0};
#ifdef _WIN32
    sprintf( cmd, "type %s && echo.", file );
#else
    sprintf( cmd, "cat %s && echo $'\n'", file );
#endif
    system( cmd );
    return 0;
}
