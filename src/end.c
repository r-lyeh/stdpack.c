#ifdef STDARC_C
#pragma once

#include <stdio.h>
#ifdef _MSC_VER
#  define ftello64 _ftelli64
#elif !defined __GNUC__
#  define ftello64 ftell
#endif

#include <stdlib.h>
#ifndef REALLOC
#define REALLOC realloc
#endif

#include <stdint.h>
#include <time.h>

static struct compressor {
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
    { LZW3, 'w', "lzw3", "lzrw3a",  lzrw3a_bounds, lzrw3a_encode, lzrw3a_decode },
    { LZSS, 's', "lzss", "lzss",    lzss_bounds, lzss_encode, lzss_decode },
    { BCM,  'B', "bcm",  "bcm",     bcm_bounds, bcm_encode, bcm_decode },
};

char *arc_nameof(unsigned flags) {
    static __thread char buf[16];
    snprintf(buf, 16, "%4s.%c", list[(flags>>4)&0x0F].name4, "0123456789ABCDEF"[flags&0xF]);
    return buf;
}

unsigned mem_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned compressor) {
    *(uint8_t*)out = compressor & 0xff;
    unsigned ret = list[(compressor >> 4) % NUM_COMPRESSORS].encode(in, inlen, (uint8_t*)out+1, outlen-1, compressor & 0x0F);
    return ret ? ret+1 : 0;
}
unsigned mem_decode(const void *in, unsigned inlen, void *out, unsigned outlen) {
    unsigned compressor = *(uint8_t*)in;
    return list[(compressor >> 4) % NUM_COMPRESSORS].decode((uint8_t*)in+1, inlen-1, out, outlen);
}
unsigned mem_bounds(unsigned inlen, unsigned compressor) {
    return list[(compressor >> 4) % NUM_COMPRESSORS].bounds(inlen, compressor & 0x0F);
}

// ---
// file options

static uint8_t STDARC_FILE_BLOCK_SIZE =  23; // 2<<(BS+12) = { 8K..256M }
static uint8_t STDARC_FILE_BLOCK_EXCESS = 0; // 16<<BE = 16, 256, 4K, 64K (16 for ulz, 256 for lpz1)

// xx yy zzzz   : 8 bits
// xx           : reserved (default = 0x11)
//    yy        : block excess [00..03] = 16<<X     = { 16, 256, 4K, 64K }
//       zzzz   : block size   [00..15] = 2<<(X+13) = { 8K..256M }

unsigned file_encode(FILE* in, FILE* out, FILE *logfile, unsigned cnum, unsigned *clist) { // multi encoder
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
    uint8_t *inbuf, *outbuf[2];
    inbuf=(uint8_t*)REALLOC(0, BS_BYTES+BE_BYTES);
    outbuf[0]=(uint8_t*)REALLOC(0, BS_BYTES*1.1+BE_BYTES);
    outbuf[1]=(uint8_t*)(cnum > 1 ? REALLOC(0, BS_BYTES*1.1+BE_BYTES) : 0);

    enum { BLOCK_PREVIEW_CHARS = 8 };
    char best_compressors_history[BLOCK_PREVIEW_CHARS+1] = {0}, best_compressors_index = BLOCK_PREVIEW_CHARS-1;
    uint8_t best = 0;

    clock_t tm = {0};
    double enctime = 0;
    if( logfile ) tm = clock();
    {
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

                    uint8_t *swap = outbuf[0];
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
    if( logfile ) enctime = (clock() - tm) / (double)CLOCKS_PER_SEC;

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

    REALLOC( outbuf[1], 0 );
    REALLOC( outbuf[0], 0 );
    REALLOC( inbuf, 0 );
    return (unsigned)total_out;
}


unsigned file_decode(FILE* in, FILE* out, FILE *logfile) { // multi decoder
    uint8_t block8; if( fread(&block8, 1,1, in ) < 1 ) return 0; 
    uint8_t excess8; if( fread(&excess8, 1,1, in ) < 1 ) return 0; 
    uint64_t BLOCK_SIZE = 1ull << block8;
    uint64_t EXCESS = 1ull << excess8;

    unsigned total = 0, outlen;
    uint8_t* inbuf=(uint8_t*)REALLOC(0, BLOCK_SIZE+EXCESS);
    uint8_t* outbuf=(uint8_t*)REALLOC(0, BLOCK_SIZE+EXCESS);

    clock_t tm = {0};
    double dectime = 0;
    if(logfile) tm = clock();
    {
        for(uint32_t inlen=0, loop=0;fread(&inlen, 1, sizeof(inlen), in)>0;++loop) {
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
            if( logfile ) fprintf(logfile, "%c\b", "\\|/-"[loop&3] );
        }
    }
    if( logfile ) dectime = (clock() - tm) / (double)CLOCKS_PER_SEC;
    if( logfile ) fprintf(logfile, "d:%.*fs ", dectime > 99 ? 1 : dectime > 9 ? 2 : 3, dectime );

    pass: goto next;
    fail: total = 0;
    next:

    REALLOC( outbuf, 0 );
    REALLOC( inbuf, 0 );
    return total;
}

#endif // STDARC_C
