// stdarc.c demo
// - rlyeh, public domain
//
// @build: cl demo.c /O2 /MT /DNDEBUG /link setargv.obj

#define STDARC_C
#include "stdarc.c"

// util
#ifdef _MSC_VER
#include <omp.h>
#define benchmark(var) for(var = -omp_get_wtime(); var < 0; var += omp_get_wtime())
#else
#include <time.h>
#define benchmark(var) for(clock_t __t = (var=-1./CLOCKS_PER_SEC, clock()); var < 0; var *= __t-clock())
#endif

// file utils
size_t file_size(const char *fname) {
    FILE *fp = fopen(fname, "rb");
    if(!fp) return 0;
    fseek(fp, 0L, SEEK_END);
    size_t sz = ftell(fp);
    fclose(fp);
    return sz;
}
char *file_read(const char *fname) {
    size_t sz = file_size(fname);
    if( !sz ) return 0;
    FILE *fp = fopen(fname, "rb");
    if(!fp) return 0;
    char *data = malloc(sz+1); data[sz] = 0;
    if( fread(data, 1,sz, fp) != sz ) { free(data); fclose(fp); return 0; }
    fclose(fp);
    return data;
}
int file_compare(FILE *in1, FILE *in2) {
    fseek(in1, 0L, SEEK_SET);
    fseek(in2, 0L, SEEK_SET);
    uint64_t sz1 = ftello64(in1);
    uint64_t sz2 = ftello64(in2);
    if(sz1 != sz2) return sz2 - sz1;
    char buffer[64*1024+1] = {0};
    while( fread(buffer,1,32*1024,in1)>0 && fread(buffer+32*1024,1,32*1024,in2)>0 ) {
        int ret = memcmp(buffer, buffer+32*1024, 32*1024);
        if( ret ) return ret;
    }
    return memcmp(buffer, buffer+32*1024, 32*1024);
}

int main(int argc, const char **argv) {
    // compression settings: high-end, middle, low-end
    unsigned list_extreme[]  = { CRSH|8, LZMA|8, BALZ|1, BCM|6 };                         // enwik8 baseline: 30..60s
    unsigned list_extra[]    = { ULZ|9, LZ4X|15, CRSH|6, DEFL|9, LZMA|7, BALZ|0, BCM|1 }; // enwik8 baseline: >12s
    unsigned list_default[]  = { ULZ|5, LZ4X|14, CRSH|4, DEFL|6, };                       // enwik8 baseline: >6s
    unsigned list_fast[]     = { ULZ|4, LZ4X|14, CRSH|0, DEFL|2, };                       // enwik8 baseline: >2s
    unsigned list_fastest[]  = { ULZ|0, LZ4X|0, PPP, LZP1, LZW3  };                       // enwik8 baseline: <1s
    unsigned list_all[]      = { ULZ,PPP,LZ4X,LZP1,LZW3,DEFL,CRSH,LZSS,BALZ,BCM,LZMA };

    if( argc <= 1 ) {
        // @todo: document everything plus --ulz,--lzma,--crsh, etc

        printf("%s %s (built: %s %s)\n\n", argv[0], STDARC_VERSION, __DATE__, __TIME__);
        printf("Usage:\n\t%s [options] files\n\n", argv[0]);
        printf("Options:\n");
        printf("\t--0..15\t\tapply compression level to selected algorithms (if possible)\n");
        printf("\t--benchmark\tbenchmark algorithms for given compression level and input files\n\n");
        printf("\t--extreme\tselect extreme compression algorithms\n");
        printf("\t--extra\t\tselect extra compression algorithms\n");
        printf("\t--default\tselect default compression algorithms (default)\n");
        printf("\t--fast\t\tselect fast compression algorithms\n");
        printf("\t--fastest\tselect fastest compression algorithms\n\n");
        for( int i = 1; i < sizeof(list)/sizeof(list[0]); ++i ) {
        printf("\t--%-8s\tselect %s compressor instead\n", list[i].name4, list[i].name);
        }
        exit(1);
    }

    unsigned *my_list = list_default, my_count = sizeof(list_default)/sizeof(list_default[0]); 

    // choose compressor list
    for( int j = 1; j < argc; ++j ) {
        /**/ if( 0 == strcmp(argv[j], "--extreme") )   my_list = list_extreme, my_count = sizeof(list_extreme)/sizeof(list_extreme[0]);
        else if( 0 == strcmp(argv[j], "--extra") )     my_list = list_extra, my_count = sizeof(list_extra)/sizeof(list_extra[0]);
        else if( 0 == strcmp(argv[j], "--fast") )      my_list = list_fast, my_count = sizeof(list_fast)/sizeof(list_fast[0]);
        else if( 0 == strcmp(argv[j], "--fastest") )   my_list = list_fastest, my_count = sizeof(list_fastest)/sizeof(list_fastest[0]);
        else if( 0 == strcmp(argv[j], "--benchmark"))  my_list = list_all, my_count = sizeof(list_all)/sizeof(list_all[0]);
    }

    // override unique compressor
    for( int j = 1; j < argc; ++j ) {
        if( argv[j][0] == '-' && argv[j][1] == '-' )
        for( int i = 1; i < NUM_COMPRESSORS - 1; ++i ) {
            if( 0 == strcmp(argv[j]+2, list[i].name4) ) {
                my_count = 1;
                my_list = &list[i].enumerator;
            }
        }
    }

    // @todo: override settings
    // BS, BE

    // patch level flags --0,--1,...--15
    for( int j = 1; j < argc; ++j ) {
        if( argv[j][0] == '-' && argv[j][1] == '-' && argv[j][2] >= '0' && argv[j][2] <= '9' ) {
            unsigned level = atoi(&argv[j][2]);
            for( int i = 0; i < my_count; ++i ) {
                my_list[i] = (my_list[i] & 0xF0) | (level > 15 ? 15 : level);
            }
        }
    }

    // benchmark
    if( my_list == list_all ) {
        double enctime, dectime, alltime = 0;
        for( int i = 1; i < argc; ++i ) {
            char *in = file_read(argv[i]);
            if( in ) {
                unsigned inlen = file_size(argv[i]);
                char *dupe = (char*)malloc(inlen + 256); // 256 == LZP1 EXCESS
                unsigned outcap = inlen * 2 + 256; // mem_bounds(inlen, flags);
                char *out = (char*)malloc(outcap);
                for( int j = 0; j < my_count; ++j ) {
                    unsigned outlen, duplen, flags = my_list[j];

                    benchmark(enctime) {
                        outlen = mem_encode(in, inlen, out, outcap, flags);
                    }
                    double ratio = outlen * 100.0 / (inlen ? inlen : 1);
                    fprintf(stdout, "\r%11u -> %11u %s %5.*f%% c:%.*fs ", 
                        inlen, outlen,
                        arc_nameof(flags),
                        ratio >= 100 ? 1 : 2, ratio,
                        enctime > 99 ? 1 : enctime > 9 ? 2 : 3, enctime);

                    memset(dupe, 0, inlen);
                    benchmark(dectime) {
                        duplen = mem_decode(out, outlen, dupe, inlen);
                    }
                    fprintf(stdout, "d:%.*fs %s\r%c\n",
                        dectime > 99 ? 1 : dectime > 9 ? 2 : 3, dectime,
                        argv[i],
                        inlen == duplen && 0 == memcmp(in, dupe, inlen) ? 'Y':'N');

                    alltime += enctime + dectime;
                }
                free(dupe);
                free(out);
                free(in);
            }
        }
        printf("%.2fs\n", alltime);
        return 0;
    }

    // de/compressor
    double enctime = 0, dectime = 0;
    uint64_t encbytes = 0, decbytes = 0, files = 0, errors = 0;
    for( int i = 1; i < argc; ++i ) {
        FILE *out = tmpfile();
        if (out) {
            FILE* in = fopen(argv[i], "rb");
            if (in) {
                unsigned x, y;
                double time;
                benchmark(time) {
                    x = file_encode(in, out, stderr, my_count, my_list);
                    encbytes += x;
                    errors += !x;
                    files += !!x;
                }
                enctime += time;

#if 0 // def NO_VERIFY
                y = (unsigned)ftell(in);
                decbytes += y;
                printf("%s\r%c\n", argv[i], "NY"[!!(x*y)]);
#else
                benchmark(time) {
                    FILE *tmp = tmpfile();
                    fseek(out, 0L, SEEK_SET);
                    y = file_decode(out, tmp, stderr);
                    decbytes += y;
                    errors += !y;
                    fseek(in, 0L, SEEK_SET);
                    fseek(tmp, 0L, SEEK_SET);
                    printf("%s\r%c\n", argv[i], "YN"[!!(errors || file_compare(tmp, in))]);
                    fclose(tmp);
                }
                dectime += time;
#endif
                fclose(in);
            }
        }
        fclose(out);
    }

    if( files > 1 ) {
        char decbuf[128] = {0};
        if( dectime > 0 ) snprintf(decbuf, 128, " d:%.*fs", dectime > 99 ? 1 : dectime > 9 ? 2 : 3, dectime);
        double ratio = (encbytes * 100.0) / (decbytes ? decbytes : 1);
        printf("%11lld -> %11lld ------ %5.*f%% c:%.*fs%s %lld files\r%c\n", 
            decbytes, encbytes,
            ratio >= 100 ? 1 : 2, ratio,
            enctime > 99 ? 1 : enctime > 9 ? 2 : 3, enctime,
            decbuf, files, "YN"[!!errors]);
    }
}
