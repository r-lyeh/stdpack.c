// stdarc.c demo
// - rlyeh, public domain
//
// @build: cl demo.c /O2 /MT /DNDEBUG /link setargv.obj
// @todo: endianness
// @todo: 0(store),1..(6)..9,10..15(uber)
// @todo: expose new/del ctx (workmem)
// @todo: compressed file seeking

// current file format:
//   header : [1<<block_size:8][1<<excess:8]
//   chunk  : [len:32] [fmt:4|lvl:4] [data:X]
//
// @todo: new format
//   header : [1<<block_size:8][1<<excess:8]
//   chunk  : [len:32][fmt|lvl:8][data:X][fmt|lvl:8][crc:32]

#ifndef STDARC_H
#define STDARC_H

#include <stdio.h>

// compressor+compression flag = 1 byte
// compressor type [0..15]: high nibble
// compression level [0..15]: low hibble

enum {
    RAW  = 0,
    PPP  = (1<<4),
    ULZ  = (2<<4),
    LZ4X = (3<<4),
    CRSH = (4<<4),
    DEFL = (5<<4),
    LZP1 = (6<<4),
    LZMA = (7<<4),
    BALZ = (8<<4),
    LZW3 = (9<<4),
    LZSS = (10<<4),
    NUM_COMPRESSORS = 12
};

unsigned file_encode(FILE* in, FILE* out, unsigned compressor); // single encoder
unsigned file_decode(FILE* in, FILE* out); // single decoder

unsigned file_encode_multi(FILE* in, FILE* out, FILE *logfile, unsigned cnum, unsigned *clist); // multi encoder
unsigned file_decode_multi(FILE* in, FILE* out, FILE *logfile); // multi decoder

int file_compare(FILE *in1, FILE *in2);

#define STDARC_VERSION "v1.0.0"

#endif

// #ifdef STDARC_C
// #pragma once

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

static unsigned raw_bounds(unsigned bytes, unsigned flags) {
    return bytes;
}
static unsigned raw_encode(const void *in, unsigned inlen, void *out, unsigned outcap, unsigned flags) {
    return memcpy(out, in, inlen), inlen;
}
static unsigned raw_decode(const void *in, unsigned inlen, void *out, unsigned outcap) {
    return memcpy(out, in, inlen), inlen;
}

struct compressor {
    // id of compressor
    unsigned enumerator;
    // name of compressor
    const char name1, *name4, *name;
    // returns worst case compression estimation for selected flags
    unsigned (*bounds)(unsigned bytes, unsigned flags);
    // returns number of bytes written. 0 if error.
    unsigned (*encode)(const void *in, unsigned inlen, void *out, unsigned outcap, unsigned flags);
    // returns number of bytes written. 0 if error.
    unsigned (*decode)(const void *in, unsigned inlen, void *out, unsigned outcap);
} list[] = {
    { RAW,  '0', "raw",  "raw",     raw_bounds, raw_encode, raw_decode },
    { PPP,  'p', "ppp",  "ppp",     ppp_bounds, ppp_encode, ppp_decode },
    { ULZ,  'u', "ulz",  "ulz",     ulz_bounds, ulz_encode, ulz_decode },
    { LZ4X, '4', "lz4x", "lz4x",    lz4x_bounds, lz4x_encode, lz4x_decode },
    { CRSH, 'c', "crsh", "crush",   crush_bounds, crush_encode, crush_decode },
    { DEFL, 'd', "defl", "deflate", deflate_bounds, deflate_encode, deflate_decode },
    { LZP1, '1', "lzp1", "lzp1",    lzp1_bounds, lzp1_encode, lzp1_decode },
    { LZMA, 'm', "lzma", "lzma",    lzma_bounds, lzma_encode, lzma_decode },
    { BALZ, 'b', "balz", "balz",    balz_bounds, balz_encode, balz_decode },
    { LZW3, 'w', "lzw3", "lzrw3-a", lzrw3a_bounds, lzrw3a_encode, lzrw3a_decode },
    { LZSS, 's', "lzss", "lzss",    lzss_bounds, lzss_encode, lzss_decode },
};

// compression settings: high-end, middle, low-end
unsigned list_extreme[]  = { CRSH|8, LZMA|8, BALZ|1 };                         // enwik8 baseline: 30..60s
unsigned list_extra[]    = { ULZ|9, LZ4X|15, CRSH|6, DEFL|9, LZMA|7, BALZ|0 }; // enwik8 baseline: >12s
unsigned list_default[]  = { ULZ|5, LZ4X|14, CRSH|4, DEFL|6,  };               // enwik8 baseline: >6s
unsigned list_fast[]     = { ULZ|4, LZ4X|14, CRSH|0, DEFL|2,  };               // enwik8 baseline: >2s
unsigned list_fastest[]  = { ULZ|0, LZ4X|0, PPP, LZP1, LZW3  };                // enwik8 baseline: <1s


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

// file options

static uint8_t STDARC_FILE_BLOCK_SIZE =  23; // 2<<(BS+12) = { 8K..256M }
static uint8_t STDARC_FILE_BLOCK_EXCESS = 0; // 16<<BE = 16, 256, 4K, 64K (16 for ulz, 256 for lpz1)

// xx yy zzzz   : 8 bits
// xx           : reserved (default = 0x11)
//    yy        : block excess [00..03] = 16<<X     = { 16, 256, 4K, 64K }
//       zzzz   : block size   [00..15] = 2<<(X+13) = { 8K..256M }


// return 0 if error

unsigned file_encode_multi(FILE* in, FILE* out, FILE *logfile, unsigned cnum, unsigned *clist) { // multi encoder
#if 0
    // uint8_t MAGIC = 0x11 << 6 | ((STDARC_FILE_BLOCK_EXCESS&3) << 4) | ((STDARC_FILE_BLOCK_SIZE-12)&15);
    // EXCESS = 16ull << ((MAGIC >> 4) & 3);
    // BLSIZE =  1ull << ((MAGIC & 15) + 13);
#else
    if( fwrite(&STDARC_FILE_BLOCK_SIZE, 1,1, out) < 1) return 0;
    if( fwrite(&STDARC_FILE_BLOCK_EXCESS, 1,1, out) < 1) return 0;
    uint64_t BS_BYTES = 1ull << STDARC_FILE_BLOCK_SIZE;
    uint64_t BE_BYTES = 1ull << STDARC_FILE_BLOCK_EXCESS;
#endif

    uint64_t total_in = 0, total_out = 0;
    uint8_t* inbuf=malloc(BS_BYTES+BE_BYTES);
    uint8_t* outbuf[2]={malloc(BS_BYTES*1.1+BE_BYTES), cnum>1 ? malloc(BS_BYTES*1.1+BE_BYTES) : 0};

    enum { BLOCK_PREVIEW_CHARS = 8 };
    char best_compressors_history[BLOCK_PREVIEW_CHARS+1] = {0}, best_compressors_index = BLOCK_PREVIEW_CHARS-1;
    uint8_t best = 0;

    double enctime;
    benchmark(enctime) {
        for( uint32_t inlen; (inlen=fread(inbuf, 1, BS_BYTES, in)) > 0 ; ) {
            uint32_t outlen[2] = {0};

            best = clist[0];
            for(unsigned i = 0; i < cnum; ++i) {
                unsigned compr = clist[i] >> 4;
                unsigned flags = clist[i] & 15;

                if(logfile) fprintf(logfile, "\r%11lld -> %11lld %4s.%c %s", (int64_t)(total_in+inlen), (int64_t)outlen[0], list[compr].name4, "0123456789ABCDEF"[flags], best_compressors_history);

                outlen[!!i] = list[compr].encode(inbuf, (unsigned)inlen, outbuf[!!i], BS_BYTES, flags);
                if(!outlen[!!i]) goto fail;

                if( i && outlen[1] < outlen[0]) {
                    best = clist[i];
                    outlen[0] = outlen[1];

                    void *swap = outbuf[0];
                    outbuf[0] = outbuf[1];
                    outbuf[1] = swap;
                }

                if(logfile) fprintf(logfile, "\r%11lld -> %11lld %4s.%c %s", (int64_t)(total_in+inlen), (int64_t)outlen[0], list[compr].name4, "0123456789ABCDEF"[flags], best_compressors_history);
            }

            uint64_t final = 4 + 1 + outlen[0]; // sizeof(outlen[0]) + sizeof(compressor) + compr data
            double ratio = final * 100.0 / (inlen ? inlen : 1);
            if(!(ratio < 97 /* && ((outlen[0] - inlen) >= 64*1024) */ )) best = 0;

            unsigned compr = best >> 4;
            unsigned flags = best & 15;

            if( compr ) {
                uint8_t packer = (compr << 4) | flags; 
                // store block length + compressor + compr data
                if( fwrite(&outlen[0], 1, 4, out) != 4 ) goto fail;
                if( fwrite(&packer, 1, 1, out) != 1 ) goto fail;
                if( fwrite(outbuf[0], 1, outlen[0], out) != outlen[0] ) goto fail;
            } else {
                uint8_t packer = 0; 
                // store block length + no-compressor + raw data
                if( fwrite(&inlen, 1, 4, out) != 4 ) goto fail;
                if( fwrite(&packer, 1, 1, out) != 1 ) goto fail;
                if( fwrite(inbuf, 1, inlen, out) != inlen ) goto fail;
            }

            total_in += inlen;
            total_out += 4 + 1 + (best ? outlen[0] : inlen);

            best_compressors_index = (best_compressors_index+1) % BLOCK_PREVIEW_CHARS;
            best_compressors_history[best_compressors_index] = list[compr].name1;
            best_compressors_history[best_compressors_index+1] = 0;
        }
    }

    if( logfile ) {
        double ratio = (total_out - 4 - 1) * 100.0 / (total_in ? total_in : 1);
        fprintf(logfile, "\r%11lld -> %11lld %4s.%c %5.*f%% c:%.*fs ", 
                total_in, total_out - 4 - 1,
                list[best>>4].name4, "0123456789ABCDEF"[best&15],
                ratio >= 100 ? 1 : 2, ratio,
                enctime > 99 ? 1 : enctime > 9 ? 2 : 3, enctime);
    }

    pass: goto next;
    fail: total_out = 0;
    next:

    free( outbuf[1] );
    free( outbuf[0] );
    free( inbuf );
    return (unsigned)total_out;
}


unsigned file_decode_multi(FILE* in, FILE* out, FILE *logfile) { // multi decoder
    uint8_t block8; if( fread(&block8, 1,1, in ) < 1 ) return 0; 
    uint8_t excess8; if( fread(&excess8, 1,1, in ) < 1 ) return 0; 
    uint64_t BLOCK_SIZE = 1ull << block8;
    uint64_t EXCESS = 1ull << excess8;

    unsigned total = 0, outlen;
    uint8_t* inbuf=malloc(BLOCK_SIZE+EXCESS);
    uint8_t* outbuf=malloc(BLOCK_SIZE+EXCESS);

    double dectime;
    benchmark(dectime) {
        for(uint32_t inlen=0;fread(&inlen, 1, sizeof(inlen), in)>0;) {
            if (inlen>(BLOCK_SIZE+EXCESS)) goto fail;

            uint8_t packer;
            if( fread(&packer, 1,sizeof(packer), in) <= 0 ) goto fail;

            if(packer) {
                // read compressed
                if (fread(inbuf, 1, inlen, in)!=inlen) goto fail;

                // decompress
                uint8_t compressor = packer >> 4;
                outlen=list[compressor % NUM_COMPRESSORS].decode(inbuf, (unsigned)inlen, outbuf, BLOCK_SIZE);
                if (!outlen) goto fail;
            } else {
                // read raw
                if (fread(outbuf, 1, inlen, in)!=inlen) goto fail;
                outlen=inlen;
            }

            if (fwrite(outbuf, 1, outlen, out) != outlen) {
                perror("fwrite() failed");
                goto fail;
            }

            total += outlen;
            if( logfile ) fprintf(logfile, "%c\b", "\\|/-"[total&3] );
        }
    }

    if( logfile ) fprintf(logfile, "d:%.*fs ", dectime > 99 ? 1 : dectime > 9 ? 2 : 3, dectime );

    pass: goto next;
    fail: total = 0;
    next:

    free( outbuf );
    free( inbuf );
    return total;
}

unsigned file_encode(FILE* in, FILE* out, unsigned compressor) { // single encoder
    return file_encode_multi(in, out, stderr, 1, &compressor);
}
unsigned file_decode(FILE* in, FILE* out) { // single decoder
    return file_decode_multi(in, out, stderr);
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

// #endif // STDARC_C


int main(int argc, const char **argv) {
    if( argc <= 1 ) {
        // @todo: document everything plus --ulz,--lzma,--crsh, etc

        printf("%s %s (built: %s %s)\n\n", argv[0], STDARC_VERSION, __DATE__, __TIME__);
        printf("Usage:\n\t%s [options] files\n\n", argv[0]);
        printf("Options:\n");
        printf("\t--0..15\t\tapply compression level to selected algorithms (if possible)\n\n");
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
        /**/ if( 0 == strcmp(argv[j], "--extreme") ) my_list = list_extreme, my_count = sizeof(list_extreme)/sizeof(list_extreme[0]);
        else if( 0 == strcmp(argv[j], "--extra") )   my_list = list_extra, my_count = sizeof(list_extra)/sizeof(list_extra[0]);
        else if( 0 == strcmp(argv[j], "--fast") )    my_list = list_fast, my_count = sizeof(list_fast)/sizeof(list_fast[0]);
        else if( 0 == strcmp(argv[j], "--fastest") ) my_list = list_fastest, my_count = sizeof(list_fastest)/sizeof(list_fastest[0]);
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

    // override flags --0,--1,...--15
    for( int j = 1; j < argc; ++j ) {
        if( argv[j][0] == '-' && argv[j][1] == '-' && argv[j][2] >= '0' && argv[j][2] <= '9' ) {
            unsigned level = atoi(&argv[j][2]);
            for( int i = 0; i < my_count; ++i ) {
                my_list[i] = (my_list[i] & 0xF0) | (level > 15 ? 15 : level);
            }
        }
    }

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
                    x = file_encode_multi(in, out, stderr, my_count, my_list);
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
                    y = file_decode_multi(out, tmp, stderr);
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
