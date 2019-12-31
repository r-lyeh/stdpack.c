// stdarc.c demo
// - rlyeh, public domain
//
// @build: cl demo.c /O2 /MT /DNDEBUG /link setargv.obj
// @todo: expose new/del ctx
// @todo: new format
//   header : [1<<block_size:8][1<<excess:8]
//   chunk  : [fmt|lvl:8][len:32][data:X][crc:32][fmt|lvl:8]

#define PPP_C
#include "ppp.c"
#define ULZ_C
#include "ulz.c"
#define LZ4X_C
#include "lz4x.c"
#define CRUSH_C
#include "crush.c"
#define DEFLATE_C
#include "deflate.c"
#define LZP1_C
#include "lzp1.c"
#define LZMA_C
#include "lzma.c"
#define BALZ_C
#include "balz.c"
#define LZRW3A_C
#include "lzrw3a.c"
#define LZSS_C
#include "lzss.c"

typedef struct compressor {
    // name of compressor
    const char* name;
    // returns worst case compression estimation for selected flags
    unsigned (*bounds)(unsigned bytes, unsigned flags);
    // returns number of bytes written. 0 if error.
    unsigned (*encode)(const void *src, unsigned src_len, void *dst, unsigned dst_cap, unsigned flags);
    // returns number of bytes written. 0 if error.
    unsigned (*decode)(const void *src, unsigned src_len, void *dst, unsigned dst_cap);
} compressor;

// ---

#include <stdint.h>

#include <stdio.h>
#ifdef _MSC_VER
#  define ftello64 _ftelli64
#endif

#ifdef _MSC_VER
#include <omp.h>
#define benchmark(var) for(var = -omp_get_wtime(); var < 0; var += omp_get_wtime())
#else
#include <time.h>
#define benchmark(var) for(clock_t __t = (var = -1, clock()); var < 0; var = (clock()-__t) / (double)CLOCKS_PER_SEC)
#endif

#define EXCESS 256          // 256 because of lzp1 // 16 because of ulz 
#define BLOCK_SIZE (1<<24)  // 16 MiB

unsigned file_encode(FILE* in, FILE* out, compressor c, unsigned flags) {
    unsigned total = 0;
    uint8_t* inbuf=malloc(BLOCK_SIZE+EXCESS);
    uint8_t* outbuf=malloc(BLOCK_SIZE+EXCESS);

    for( uint32_t inlen; (inlen=fread(inbuf, 1, BLOCK_SIZE, in)) > 0 ; ) {
        const uint32_t outlen = c.encode(inbuf, (unsigned)inlen, outbuf, BLOCK_SIZE, flags);

        // @todo: store uncompressed on errors
        if( outlen <= 0 ) { total = 0; break; }

        // @todo: add checks
        fwrite(&outlen, 1, sizeof(outlen), out);
        fwrite(outbuf, 1, outlen, out);

        total += outlen;
        fprintf(stderr, "%11lld -> %11lld\r", ftello64(in), ftello64(out));
    }

    free( outbuf );
    free( inbuf );
    return total;
}

unsigned file_decode(FILE* in, FILE* out, compressor c, unsigned flags) {
    unsigned total = 0;
    uint8_t* inbuf=malloc(BLOCK_SIZE+EXCESS);
    uint8_t* outbuf=malloc(BLOCK_SIZE+EXCESS);

    for(uint32_t inlen=0;fread(&inlen, 1, sizeof(inlen), in)>0;) {
        if (inlen>(BLOCK_SIZE+EXCESS) || fread(inbuf, 1, inlen, in)!=inlen)
            { total = 0; break; }

        unsigned outlen=c.decode(inbuf, (unsigned)inlen, outbuf, BLOCK_SIZE);
        if (!outlen)
            { total = 0; break; }

        if (fwrite(outbuf, 1, outlen, out) != outlen) {
            perror("fwrite() failed");
            total = 0;
            break;
            exit(1);
        }

        total += outlen;
        fprintf(stderr, "%11lld -> %11lld\r", ftello64(in), ftello64(out));
    }

    free( outbuf );
    free( inbuf );
    return total;
}

unsigned file_tests(const char *infile, compressor cs, unsigned flags) {
    double enctime = 0, dectime = 0;

    FILE *in = fopen(infile, "rb");
    if( !in ) {
        return 0;
    }

    FILE *tmp = tmpfile(), *out = tmpfile();

        benchmark(enctime) {
            file_encode(in, tmp, cs, flags);
        }
        double zratio = ftello64(tmp) * 100.0 / (1+ftello64(in));
        int is_compressable = zratio < 97; // && ((zlen - ulen) >= 64*1024);
        printf("%11lld -> %11lld %-4s %5.2f%% c:%.*fs d:%1.3fs %s\r", 
                ftello64(in), ftello64(tmp),
                cs.name, zratio,
                enctime > 10 ? 2 : 3, enctime, dectime, infile);

        fseek(tmp, 0L, SEEK_SET);
        benchmark(dectime) {
            file_decode(tmp, out, cs, flags);
        }
        printf("%11lld -> %11lld %-4s %5.2f%% c:%.*fs d:%1.3fs %s\r", 
                ftello64(out), ftello64(tmp),
                cs.name, zratio,
                enctime > 99 ? 1 : enctime > 9 ? 2 : 3, enctime, dectime, infile);

        fseek(in, 0L, SEEK_SET);
        fseek(out, 0L, SEEK_SET);
        int errors = 0;
        char buffer[64*1024+1] = {0};
        while( fread(buffer,1,32*1024,in)>0 && fread(buffer+32*1024,1,32*1024,out)>0 ) {
            errors += !!memcmp(buffer, buffer+32*1024, 32*1024);
        }
        errors += !!memcmp(buffer, buffer+32*1024, 32*1024);
        printf("\r%s\n", errors ? "N":"Y");

    fseek(tmp, 0L, SEEK_END);
    size_t zlen = ftell(tmp);

    fclose(in);
    fclose(out);
    fclose(tmp);

    return errors ? 0 : (unsigned)zlen;
}


int main(int argc, const char **argv) {
    if( argc < 2 ) {
        return -printf("%s [--0..9] file\n", argv[0]);
    }

    compressor list[] = {
        { "ppp",  ppp_bounds, ppp_encode, ppp_decode },
        { "lzp1", lzp1_bounds, lzp1_encode, lzp1_decode },
        { "lzss", lzss_bounds, lzss_encode, lzss_decode },
        { "lzrw", lzrw3a_bounds, lzrw3a_encode, lzrw3a_decode },
        { "lz4x", lz4x_bounds, lz4x_encode, lz4x_decode },
        { "ulz",  ulz_bounds, ulz_encode, ulz_decode },
        { "defl", deflate_bounds, deflate_encode, deflate_decode },
        { "crsh", crush_bounds, crush_encode, crush_decode },
        { "balz", balz_bounds, balz_encode, balz_decode },
        { "lzma", lzma_bounds, lzma_encode, lzma_decode },
    };

    unsigned level = argv[1][0]=='-' && argv[1][1]=='-' ? atoi(&argv[1][2]) : 7;

    double elapsed;
    benchmark(elapsed) {
        for( int j = 0; j < sizeof(list) / sizeof(list[0]); ++j ) {
            int zlen = 0, num_files = 0;
            for( int i = 1; i < argc; ++i ) {
                int rc = file_tests(argv[i], list[j], level);
                num_files += rc > 0;
                zlen += rc > 0 ? rc : 0;
            }
            if( num_files > 1 )
            printf("%4d files, %11d bytes (%s)\n", num_files, zlen, list[j].name);
        }
    }
    printf("%.2fs (level=%d)\n", elapsed, level);
}
