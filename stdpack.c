#line 1 "amalgamated_pack.h" 
// stdpack.c
// - rlyeh, public domain
//
// current file format:
//   header : [1<<block_size:8][1<<excess:8]
//   chunk  : [len:32] [fmt:4|lvl:4] [data:X]
//
// @todo: new format
//   header : [1<<block_size:8][1<<excess:8]
//   chunk  : [len:32][fmt|lvl:8][data:X][fmt|lvl:8][crc:32]
//
// @todo: endianness
// @todo: 0(store),1..(6)..9,10..15(uber)
// @todo: expose new/del ctx (workmem)
// @todo: compressed file seeking

#ifndef STDPACK_H
#define STDPACK_H
#define STDPACK_VERSION "v1.0.0"

#include <stdio.h>

// compressor type [0..15]: high nibble
// compression level/flags [0..15]: low hibble
// compressor_type << 4 + compression_level = 1 byte

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
    BCM  = (11<<4),
    NUM_COMPRESSORS = 13
};

// mem de/encoder
unsigned mem_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned compressor);
unsigned mem_decode(const void *in, unsigned inlen, void *out, unsigned outlen);
unsigned mem_bounds(unsigned inlen, unsigned compressor);

// file de/encoder
unsigned file_encode(FILE* in, FILE* out, FILE *logfile, unsigned cnum, unsigned *clist);
unsigned file_decode(FILE* in, FILE* out, FILE *logfile);

#endif

#ifdef STDPACK_C
#pragma once
#define RAW_C
#define PPP_C
#define ULZ_C
#define LZ4X_C
#define CRUSH_C
#define DEFLATE_C
#define LZP1_C
#define LZMA_C
#define BALZ_C
#define LZRW3A_C
#define LZSS_C
#define BCM_C
#define ZIP_C
#define TAR_C
#define PAK_C
#define VFS_C
#define DIR_C
#endif

#line 1 "amalgamated_balz.c" 
// balz.cpp is written and placed in the public domain by Ilya Muravyov
// additional code by @r-lyeh (public domain)

unsigned balz_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags /*[0..1]*/);
unsigned balz_decode(const void *in, unsigned inlen, void *out, unsigned outlen);
unsigned balz_bounds(unsigned inlen, unsigned flags);


#ifdef BALZ_C
#pragma once
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_DISABLE_PERFCRIT_LOCKS
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct mfile {
    uint8_t *begin, *seek, *end;
} mfile;
int minit(mfile *f, const void *ptr, int len) {
    f->begin = f->seek = f->end = (uint8_t*)ptr;
    f->end += len;
    return 0;
}
int mread(mfile *m, void *buf, int len) {
    if( len >= (m->end - m->seek) ) len = (m->end - m->seek);
    memcpy(buf,m->seek,len); m->seek += len;
    return len;
}
int mwrite(mfile *m, const void *buf, int len) {
    if( len >= (m->end - m->seek) ) len = (m->end - m->seek);
    memcpy(m->seek,buf,len); m->seek += len;
    return len;
}
int mtell(mfile *m) {
    return m->seek - m->begin;
}
int mavail(mfile *m) {
    return m->end - m->seek;
}
int mputc(mfile *m, int i) {
    uint8_t ch = i;
    return mwrite(m, &ch, 1);
}
int mgetc(mfile *m) { 
    if( mavail(m) <= 0 ) return -1;
    uint8_t ch; mread(m, &ch, 1); return ch;
}




typedef struct Counter {
    uint16_t p1;
    uint16_t p2;
} Counter;

void CounterCtor(Counter *c) {
    c->p1 = 1<<15;
    c->p2 = 1<<15;
}

uint32_t CounterP(const Counter *c) {
    return c->p1+c->p2;
}

void CounterUpdate0(Counter *c) {
    c->p1-=c->p1>>3;
    c->p2-=c->p2>>6;
}

void CounterUpdate1(Counter *c) {
    c->p1+=(c->p1^65535)>>3;
    c->p2+=(c->p2^65535)>>6;
}

typedef struct Encoder {
    uint32_t code;
    uint32_t low;
    uint32_t high;
    mfile *in, *out;
} Encoder;

void EncoderCtor(Encoder *e, mfile *in, mfile *out) {
    e->code = e->low = 0; e->high = -1;
    e->in = in;
    e->out = out;
}

void EncoderEncode(Encoder *e, int bit, Counter *c) {
    const uint32_t mid=e->low+((((uint64_t)e->high-e->low)*(CounterP(c)<<15))>>32);

    if (bit) {
        e->high=mid;
        CounterUpdate1(c);
    } else {
        e->low=mid+1;
        CounterUpdate0(c);
    }

    while ((e->low^e->high)<(1<<24)) {
        mputc(e->out, e->low>>24);
        e->low<<=8;
        e->high=(e->high<<8)|255;
    }
}

void EncoderFlush(Encoder *e) {
    for (int i=0; i<4; ++i) {
        mputc(e->out, e->low>>24);
        e->low<<=8;
    }
}

void EncoderInit(Encoder *e) {
    for (int i=0; i<4; ++i)
        e->code=(e->code<<8)|mgetc(e->in);
}

int EncoderDecode(Encoder *e, Counter *c) {
    const uint32_t mid=e->low+((((uint64_t)e->high-e->low)*(CounterP(c)<<15))>>32);

    const int bit=(e->code<=mid);
    if (bit) {
        e->high=mid;
        CounterUpdate1(c);
    } else {
        e->low=mid+1;
        CounterUpdate0(c);
    }

    while ((e->low^e->high)<(1<<24)) {
        e->code=(e->code<<8)|mgetc(e->in);
        e->low<<=8;
        e->high=(e->high<<8)|255;
    }

    return bit;
}

enum { BALZ_TAB_BITS=7 };
enum { BALZ_TAB_SIZE=1<<BALZ_TAB_BITS };
enum { BALZ_TAB_MASK=BALZ_TAB_SIZE-1 };

typedef struct CM {
    Encoder encoder;
    Counter counter1[256][512];
    Counter counter2[256][BALZ_TAB_SIZE];
} CM;

void CMCtor(CM *cm, mfile *in, mfile *out) {
    EncoderCtor(&cm->encoder, in, out);
    for( int i = 0; i < 256; ++i) for( int j = 0; j < 512;           ++j) CounterCtor(&cm->counter1[i][j]);
    for( int i = 0; i < 256; ++i) for( int j = 0; j < BALZ_TAB_SIZE; ++j) CounterCtor(&cm->counter2[i][j]);
}

void CMInit(CM *cm) {
    EncoderInit(&cm->encoder);
}

void CMEncode(CM *cm, int t, int c1) {
    int ctx=1;
    while (ctx<512) {
        const int bit=((t&256)!=0);
        t+=t;
        EncoderEncode(&cm->encoder, bit, &cm->counter1[c1][ctx]);
        ctx+=ctx+bit;
    }
}

void CMEncodeIdx(CM *cm, int x, int c2) {
    int ctx=1;
    while (ctx<BALZ_TAB_SIZE) {
        const int bit=((x&(BALZ_TAB_SIZE>>1))!=0);
        x+=x;
        EncoderEncode(&cm->encoder, bit, &cm->counter2[c2][ctx]);
        ctx+=ctx+bit;
    }
}

int CMDecode(CM *cm, int c1) {
    int ctx=1;
    while (ctx<512)
        ctx+=ctx+EncoderDecode(&cm->encoder, &cm->counter1[c1][ctx]);
    return ctx-512;
}

int CMDecodeIdx(CM *cm, int c2) {
    int ctx=1;
    while (ctx<BALZ_TAB_SIZE)
        ctx+=ctx+EncoderDecode(&cm->encoder, &cm->counter2[c2][ctx]);
    return ctx-BALZ_TAB_SIZE;
}

enum { BALZ_MIN_MATCH=3 };
enum { BALZ_MAX_MATCH=255+BALZ_MIN_MATCH };

enum { BALZ_BUF_BITS=25 };
enum { BALZ_BUF_SIZE=1<<BALZ_BUF_BITS };
enum { BALZ_BUF_MASK=BALZ_BUF_SIZE-1 };

uint8_t buf[BALZ_BUF_SIZE];
uint32_t tab[1<<16][BALZ_TAB_SIZE];
int cnt[1<<16];

void balz_init() {
    size_t buflen = sizeof(uint8_t)  * BALZ_BUF_SIZE;
    memset(buf, 0, buflen);

    size_t cntlen = sizeof(int)      * (1<<16);
    memset(cnt, 0, cntlen);

    for(int i = 0; i < (1<<16); ++i)
    for(int j = 0; j < BALZ_TAB_SIZE; ++j)
    tab[i][j] = 0;
}

// E8E9 preprocessor to improve compression of x86 (EXE and DLL) files.
// The preprocessor replaces relative CALL and JMP addresses with absolute addresses,
// which improves compression because an address may appear multiple times.
// Many other compressors use this technique.

#define e8e9_transform(FWD,n) do { \
    const int end=n-8; \
    int p=0; \
    while((*((int*)&buf[p])!=0x4550)&&(++p<end)); /* unaligned */ \
    while (p<end) { \
        if ((buf[p++]&254)==0xe8) { \
            int *addr=(int*)&buf[p]; /* unaligned */ \
            if (FWD) { \
                if ((*addr>=-p)&&(*addr<(n-p))) \
                    *addr+=p; \
                else if ((*addr>0)&&(*addr<n)) \
                    *addr-=n; \
            } else { \
                if (*addr<0) { \
                    if ((*addr+p)>=0) \
                        *addr+=n; \
                } \
                else if (*addr<n) \
                    *addr-=p; \
            } \
            p+=4; \
        } \
    } \
} while(0)

static inline uint32_t get_hash(int p) {
    return (((*(uint32_t*)(&buf[p]))&0xffffff) *2654435769UL)&~BALZ_BUF_MASK; // Little-endian+unaligned
}

static inline int get_pts(int len, int x) {
    return len>=BALZ_MIN_MATCH?(len<<BALZ_TAB_BITS)-x:((BALZ_MIN_MATCH-1)<<BALZ_TAB_BITS)-8;
}

int get_pts_at(int p, int n) {
    const int c2=*(uint16_t*)&buf[p-2]; // unaligned
    const uint32_t hash=get_hash(p);

    int len=BALZ_MIN_MATCH-1;
    int idx=BALZ_TAB_SIZE;

    int max_match=n-p;
    if (max_match>BALZ_MAX_MATCH)
        max_match=BALZ_MAX_MATCH;

    for (int x=0; x<BALZ_TAB_SIZE; ++x) {
        const uint32_t d=tab[c2][(cnt[c2]-x)&BALZ_TAB_MASK];
        if (!d)
            break;

        if ((d&~BALZ_BUF_MASK)!=hash)
            continue;

        const int s=d&BALZ_BUF_MASK;
        if ((buf[s+len]!=buf[p+len])||(buf[s]!=buf[p]))
            continue;

        int l=0;
        while (++l<max_match)
            if (buf[s+l]!=buf[p+l])
                break;
        if (l>len) {
            idx=x;
            len=l;
            if (l==max_match)
                break;
        }
    }

    return get_pts(len, idx);
}

int balz_compress(const uint8_t *in, unsigned inlen, uint8_t *out, unsigned outlen, unsigned is_max) {
    balz_init();

    *out++ = (inlen >> 24) & 255;
    *out++ = (inlen >> 16) & 255;
    *out++ = (inlen >>  8) & 255;
    *out++ = (inlen >>  0) & 255;
    outlen -= 4;

    mfile inf, outf;
    minit(&inf, in, inlen);
    minit(&outf, out, outlen);

    CM cm;
    CMCtor(&cm, &inf, &outf);

    int best_idx[BALZ_MAX_MATCH+1];

    int n;
    while ((n=mread(&inf, buf, BALZ_BUF_SIZE))>0) {
        //e8e9_transform(1,n);

        memset(tab, 0, sizeof(tab));

        int p=0;

        while ((p<2)&&(p<n))
            CMEncode(&cm, buf[p++], 0);

        while (p<n) {
            const int c2=*(uint16_t*)&buf[p-2]; // unaligned
            const uint32_t hash=get_hash(p);

            int len=BALZ_MIN_MATCH-1;
            int idx=BALZ_TAB_SIZE;

            int max_match=n-p;
            if (max_match>BALZ_MAX_MATCH)
                max_match=BALZ_MAX_MATCH;

            for (int x=0; x<BALZ_TAB_SIZE; ++x) {
                const uint32_t d=tab[c2][(cnt[c2]-x)&BALZ_TAB_MASK];
                if (!d)
                    break;

                if ((d&~BALZ_BUF_MASK)!=hash)
                    continue;

                const int s=d&BALZ_BUF_MASK;
                if ((buf[s+len]!=buf[p+len])||(buf[s]!=buf[p]))
                    continue;

                int l=0;
                while (++l<max_match)
                    if (buf[s+l]!=buf[p+l])
                        break;
                if (l>len) {
                    for (int i=l; i>len; --i)
                        best_idx[i]=x;
                    idx=x;
                    len=l;
                    if (l==max_match)
                        break;
                }
            }

            if ((is_max)&&(len>=BALZ_MIN_MATCH)) {
                int sum=get_pts(len, idx)+get_pts_at(p+len, n);

                if (sum<get_pts(len+BALZ_MAX_MATCH, 0)) {
                    const int lookahead=len;

                    for (int i=1; i<lookahead; ++i) {
                        const int tmp=get_pts(i, best_idx[i])+get_pts_at(p+i, n);
                        if (tmp>sum) {
                            sum=tmp;
                            len=i;
                        }
                    }

                    idx=best_idx[len];
                }
            }

            tab[c2][++cnt[c2]&BALZ_TAB_MASK]=hash|p;

            if (len>=BALZ_MIN_MATCH) {
                CMEncode(&cm, (256-BALZ_MIN_MATCH)+len, buf[p-1]);
                CMEncodeIdx(&cm, idx, buf[p-2]);
                p+=len;
            } else {
                CMEncode(&cm, buf[p], buf[p-1]);
                ++p;
            }
        }
    }

    EncoderFlush(&cm.encoder);

    if ( (inf.seek - inf.begin) != inlen) {
        return 0; // size mismatch error
    }

    return (int)(outf.seek - outf.begin) + 4;
}

int balz_decompress(const uint8_t *in, unsigned inlen, uint8_t *out, unsigned outlen) {
    balz_init();

    uint32_t flen32 = 0;
    flen32 |= ((uint32_t)*in++) << 24;
    flen32 |= ((uint32_t)*in++) << 16;
    flen32 |= ((uint32_t)*in++) <<  8;
    flen32 |= ((uint32_t)*in++) <<  0;
    outlen = flen32;
    int flen = flen32; inlen -= 4;

    mfile inf, outf;
    minit(&inf, in, inlen);
    minit(&outf, out, outlen);

    CM cm;
    CMCtor(&cm, &inf, &outf);
    CMInit(&cm);

    #define balz_src_avail ((int)(inf.end - inf.seek))
    #define balz_dst_avail ((int)(outf.end - outf.seek))
    #define balz_dst_written ((int)(outf.seek - outf.begin))

    while(/*(balz_src_avail > 0) &&*/ (balz_dst_written != flen)) {
        int p=0;

        while ((p<2) && ((p+balz_dst_written)<flen)) {
            const int t=CMDecode(&cm, 0);
            if (t>=256) {
                return 0; // corrupt file error
            }
            buf[p++]=t;
        }

        while ((p < BALZ_BUF_SIZE) && (p+balz_dst_written<flen)) { // (balz_src_avail > 0)) {
            const int tmp=p;
            const int c2=*(uint16_t*)(&buf[p-2]); // unaligned

            const int t=CMDecode(&cm, buf[p-1]);
            if (t>=256) {
                int len=t-256;
                int s=tab[c2][(cnt[c2]-CMDecodeIdx(&cm, buf[p-2]))&BALZ_TAB_MASK];

                buf[p++]=buf[s++];
                buf[p++]=buf[s++];
                buf[p++]=buf[s++];
                while (len--)
                    buf[p++]=buf[s++];
            }
            else
                buf[p++]=t;

            tab[c2][++cnt[c2]&BALZ_TAB_MASK]=tmp;
        }

        //e8e9_transform(0,p);

        mwrite(&outf, buf, p);
    }

    return (int)(outf.seek - outf.begin);
}

unsigned balz_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags /*[0..1]*/) {
    unsigned level = flags > 0 ? 1 : 0;
    return (unsigned)balz_compress((const uint8_t *)in, inlen, (uint8_t*)out, outlen, level);
}
unsigned balz_decode(const void *in, unsigned inlen, void *out, unsigned outlen) {
    return (unsigned)balz_decompress((const uint8_t *)in, inlen, (uint8_t*)out, outlen);
}
unsigned balz_bounds(unsigned inlen, unsigned flags) {
    return (unsigned)(inlen * 1.1) + 16; // @todo: check src
}

#endif // BALZ_C

#ifdef BALZ_DEMO
#pragma once
#include <string.h>
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level=1;
    char out[128];
    unsigned outlen = balz_encode(longcopy, strlen(longcopy)+1, out, 128, level);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    unsigned unpacked = balz_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // BALZ_DEMO

#line 1 "amalgamated_bcm_bwt.c" 
#ifndef BCM_C
// do nothing
#elif defined BCM_NO_ENCODER
// dummy
int bcm_divbwt(const unsigned char *T, unsigned char *U, int *A, int n) { return -1; }
#else
/*
 * divsufsort.h for libdivsufsort-lite
 * Copyright (c) 2003-2008 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _DIVSUFSORT_H
#define _DIVSUFSORT_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*- Prototypes -*/

/**
 * Constructs the suffix array of a given string.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-1] The output array of suffixes.
 * @param n The length of the given string.
 * @return 0 if no error occurred, -1 or -2 otherwise.
 */
int
bcm_divsufsort(const unsigned char *T, int *SA, int n);

/**
 * Constructs the burrows-wheeler transformed string of a given string.
 * @param T[0..n-1] The input string.
 * @param U[0..n-1] The output string. (can be T)
 * @param A[0..n-1] The temporary array. (can be NULL)
 * @param n The length of the given string.
 * @return The primary index if no error occurred, -1 or -2 otherwise.
 */
int
bcm_divbwt(const unsigned char *T, unsigned char *U, int *A, int n);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _DIVSUFSORT_H */

/*
 * divsufsort.c for libdivsufsort-lite
 * Copyright (c) 2003-2008 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
# include <omp.h>
#endif
//#include "bcm_divsufsort.h"


/*- Constants -*/
#define INLINE __inline
#if defined(ALPHABET_SIZE) && (ALPHABET_SIZE < 1)
# undef ALPHABET_SIZE
#endif
#if !defined(ALPHABET_SIZE)
# define ALPHABET_SIZE (256)
#endif
#define BUCKET_A_SIZE (ALPHABET_SIZE)
#define BUCKET_B_SIZE (ALPHABET_SIZE * ALPHABET_SIZE)
#if defined(SS_INSERTIONSORT_THRESHOLD)
# if SS_INSERTIONSORT_THRESHOLD < 1
#  undef SS_INSERTIONSORT_THRESHOLD
#  define SS_INSERTIONSORT_THRESHOLD (1)
# endif
#else
# define SS_INSERTIONSORT_THRESHOLD (8)
#endif
#if defined(SS_BLOCKSIZE)
# if SS_BLOCKSIZE < 0
#  undef SS_BLOCKSIZE
#  define SS_BLOCKSIZE (0)
# elif 32768 <= SS_BLOCKSIZE
#  undef SS_BLOCKSIZE
#  define SS_BLOCKSIZE (32767)
# endif
#else
# define SS_BLOCKSIZE (1024)
#endif
/* minstacksize = log(SS_BLOCKSIZE) / log(3) * 2 */
#if SS_BLOCKSIZE == 0
# define SS_MISORT_STACKSIZE (96)
#elif SS_BLOCKSIZE <= 4096
# define SS_MISORT_STACKSIZE (16)
#else
# define SS_MISORT_STACKSIZE (24)
#endif
#define SS_SMERGE_STACKSIZE (32)
#define TR_INSERTIONSORT_THRESHOLD (8)
#define TR_STACKSIZE (64)


/*- Macros -*/
#ifndef SWAP
# define SWAP(_a, _b) do { t = (_a); (_a) = (_b); (_b) = t; } while(0)
#endif /* SWAP */
#ifndef MIN
# define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif /* MIN */
#ifndef MAX
# define MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))
#endif /* MAX */
#define STACK_PUSH(_a, _b, _c, _d)\
    do {\
        assert(ssize < STACK_SIZE);\
        stack[ssize].a = (_a), stack[ssize].b = (_b),\
        stack[ssize].c = (_c), stack[ssize++].d = (_d);\
    } while(0)
#define STACK_PUSH5(_a, _b, _c, _d, _e)\
    do {\
        assert(ssize < STACK_SIZE);\
        stack[ssize].a = (_a), stack[ssize].b = (_b),\
        stack[ssize].c = (_c), stack[ssize].d = (_d), stack[ssize++].e = (_e);\
    } while(0)
#define STACK_POP(_a, _b, _c, _d)\
    do {\
        assert(0 <= ssize);\
        if(ssize == 0) { return; }\
        (_a) = stack[--ssize].a, (_b) = stack[ssize].b,\
        (_c) = stack[ssize].c, (_d) = stack[ssize].d;\
    } while(0)
#define STACK_POP5(_a, _b, _c, _d, _e)\
    do {\
        assert(0 <= ssize);\
        if(ssize == 0) { return; }\
        (_a) = stack[--ssize].a, (_b) = stack[ssize].b,\
        (_c) = stack[ssize].c, (_d) = stack[ssize].d, (_e) = stack[ssize].e;\
    } while(0)
#define BUCKET_A(_c0) bucket_A[(_c0)]
#if ALPHABET_SIZE == 256
#define BUCKET_B(_c0, _c1) (bucket_B[((_c1) << 8) | (_c0)])
#define BUCKET_BSTAR(_c0, _c1) (bucket_B[((_c0) << 8) | (_c1)])
#else
#define BUCKET_B(_c0, _c1) (bucket_B[(_c1) * ALPHABET_SIZE + (_c0)])
#define BUCKET_BSTAR(_c0, _c1) (bucket_B[(_c0) * ALPHABET_SIZE + (_c1)])
#endif


/*- Private Functions -*/

static const int lg_table[256]= {
 -1,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7
};

#if (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE)

static INLINE
int
ss_ilg(int n) {
#if SS_BLOCKSIZE == 0
    return (n & 0xffff0000) ?
                    ((n & 0xff000000) ?
                        24 + lg_table[(n >> 24) & 0xff] :
                        16 + lg_table[(n >> 16) & 0xff]) :
                    ((n & 0x0000ff00) ?
                         8 + lg_table[(n >>  8) & 0xff] :
                         0 + lg_table[(n >>  0) & 0xff]);
#elif SS_BLOCKSIZE < 256
    return lg_table[n];
#else
    return (n & 0xff00) ?
                    8 + lg_table[(n >> 8) & 0xff] :
                    0 + lg_table[(n >> 0) & 0xff];
#endif
}

#endif /* (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE) */

#if SS_BLOCKSIZE != 0

static const int sqq_table[256] = {
    0,  16,  22,  27,  32,  35,  39,  42,  45,  48,  50,  53,  55,  57,  59,  61,
 64,  65,  67,  69,  71,  73,  75,  76,  78,  80,  81,  83,  84,  86,  87,  89,
 90,  91,  93,  94,  96,  97,  98,  99, 101, 102, 103, 104, 106, 107, 108, 109,
110, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
128, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
143, 144, 144, 145, 146, 147, 148, 149, 150, 150, 151, 152, 153, 154, 155, 155,
156, 157, 158, 159, 160, 160, 161, 162, 163, 163, 164, 165, 166, 167, 167, 168,
169, 170, 170, 171, 172, 173, 173, 174, 175, 176, 176, 177, 178, 178, 179, 180,
181, 181, 182, 183, 183, 184, 185, 185, 186, 187, 187, 188, 189, 189, 190, 191,
192, 192, 193, 193, 194, 195, 195, 196, 197, 197, 198, 199, 199, 200, 201, 201,
202, 203, 203, 204, 204, 205, 206, 206, 207, 208, 208, 209, 209, 210, 211, 211,
212, 212, 213, 214, 214, 215, 215, 216, 217, 217, 218, 218, 219, 219, 220, 221,
221, 222, 222, 223, 224, 224, 225, 225, 226, 226, 227, 227, 228, 229, 229, 230,
230, 231, 231, 232, 232, 233, 234, 234, 235, 235, 236, 236, 237, 237, 238, 238,
239, 240, 240, 241, 241, 242, 242, 243, 243, 244, 244, 245, 245, 246, 246, 247,
247, 248, 248, 249, 249, 250, 250, 251, 251, 252, 252, 253, 253, 254, 254, 255
};

static INLINE
int
ss_isqrt(int x) {
    int y, e;

    if(x >= (SS_BLOCKSIZE * SS_BLOCKSIZE)) { return SS_BLOCKSIZE; }
    e = (x & 0xffff0000) ?
                ((x & 0xff000000) ?
                    24 + lg_table[(x >> 24) & 0xff] :
                    16 + lg_table[(x >> 16) & 0xff]) :
                ((x & 0x0000ff00) ?
                     8 + lg_table[(x >>  8) & 0xff] :
                     0 + lg_table[(x >>  0) & 0xff]);

    if(e >= 16) {
        y = sqq_table[x >> ((e - 6) - (e & 1))] << ((e >> 1) - 7);
        if(e >= 24) { y = (y + 1 + x / y) >> 1; }
        y = (y + 1 + x / y) >> 1;
    } else if(e >= 8) {
        y = (sqq_table[x >> ((e - 6) - (e & 1))] >> (7 - (e >> 1))) + 1;
    } else {
        return sqq_table[x] >> 4;
    }

    return (x < (y * y)) ? y - 1 : y;
}

#endif /* SS_BLOCKSIZE != 0 */


/*---------------------------------------------------------------------------*/

/* Compares two suffixes. */
static INLINE
int
ss_compare(const unsigned char *T,
                     const int *p1, const int *p2,
                     int depth) {
    const unsigned char *U1, *U2, *U1n, *U2n;

    for(U1 = T + depth + *p1,
            U2 = T + depth + *p2,
            U1n = T + *(p1 + 1) + 2,
            U2n = T + *(p2 + 1) + 2;
            (U1 < U1n) && (U2 < U2n) && (*U1 == *U2);
            ++U1, ++U2) {
    }

    return U1 < U1n ?
                (U2 < U2n ? *U1 - *U2 : 1) :
                (U2 < U2n ? -1 : 0);
}


/*---------------------------------------------------------------------------*/

#if (SS_BLOCKSIZE != 1) && (SS_INSERTIONSORT_THRESHOLD != 1)

/* Insertionsort for small size groups */
static
void
ss_insertionsort(const unsigned char *T, const int *PA,
                                 int *first, int *last, int depth) {
    int *i, *j;
    int t;
    int r;

    for(i = last - 2; first <= i; --i) {
        for(t = *i, j = i + 1; 0 < (r = ss_compare(T, PA + t, PA + *j, depth));) {
            do { *(j - 1) = *j; } while((++j < last) && (*j < 0));
            if(last <= j) { break; }
        }
        if(r == 0) { *j = ~*j; }
        *(j - 1) = t;
    }
}

#endif /* (SS_BLOCKSIZE != 1) && (SS_INSERTIONSORT_THRESHOLD != 1) */


/*---------------------------------------------------------------------------*/

#if (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE)

static INLINE
void
ss_fixdown(const unsigned char *Td, const int *PA,
                     int *SA, int i, int size) {
    int j, k;
    int v;
    int c, d, e;

    for(v = SA[i], c = Td[PA[v]]; (j = 2 * i + 1) < size; SA[i] = SA[k], i = k) {
        d = Td[PA[SA[k = j++]]];
        if(d < (e = Td[PA[SA[j]]])) { k = j; d = e; }
        if(d <= c) { break; }
    }
    SA[i] = v;
}

/* Simple top-down heapsort. */
static
void
ss_heapsort(const unsigned char *Td, const int *PA, int *SA, int size) {
    int i, m;
    int t;

    m = size;
    if((size % 2) == 0) {
        m--;
        if(Td[PA[SA[m / 2]]] < Td[PA[SA[m]]]) { SWAP(SA[m], SA[m / 2]); }
    }

    for(i = m / 2 - 1; 0 <= i; --i) { ss_fixdown(Td, PA, SA, i, m); }
    if((size % 2) == 0) { SWAP(SA[0], SA[m]); ss_fixdown(Td, PA, SA, 0, m); }
    for(i = m - 1; 0 < i; --i) {
        t = SA[0], SA[0] = SA[i];
        ss_fixdown(Td, PA, SA, 0, i);
        SA[i] = t;
    }
}


/*---------------------------------------------------------------------------*/

/* Returns the median of three elements. */
static INLINE
int *
ss_median3(const unsigned char *Td, const int *PA,
                     int *v1, int *v2, int *v3) {
    int *t;
    if(Td[PA[*v1]] > Td[PA[*v2]]) { SWAP(v1, v2); }
    if(Td[PA[*v2]] > Td[PA[*v3]]) {
        if(Td[PA[*v1]] > Td[PA[*v3]]) { return v1; }
        else { return v3; }
    }
    return v2;
}

/* Returns the median of five elements. */
static INLINE
int *
ss_median5(const unsigned char *Td, const int *PA,
                     int *v1, int *v2, int *v3, int *v4, int *v5) {
    int *t;
    if(Td[PA[*v2]] > Td[PA[*v3]]) { SWAP(v2, v3); }
    if(Td[PA[*v4]] > Td[PA[*v5]]) { SWAP(v4, v5); }
    if(Td[PA[*v2]] > Td[PA[*v4]]) { SWAP(v2, v4); SWAP(v3, v5); }
    if(Td[PA[*v1]] > Td[PA[*v3]]) { SWAP(v1, v3); }
    if(Td[PA[*v1]] > Td[PA[*v4]]) { SWAP(v1, v4); SWAP(v3, v5); }
    if(Td[PA[*v3]] > Td[PA[*v4]]) { return v4; }
    return v3;
}

/* Returns the pivot element. */
static INLINE
int *
ss_pivot(const unsigned char *Td, const int *PA, int *first, int *last) {
    int *middle;
    int t;

    t = last - first;
    middle = first + t / 2;

    if(t <= 512) {
        if(t <= 32) {
            return ss_median3(Td, PA, first, middle, last - 1);
        } else {
            t >>= 2;
            return ss_median5(Td, PA, first, first + t, middle, last - 1 - t, last - 1);
        }
    }
    t >>= 3;
    first  = ss_median3(Td, PA, first, first + t, first + (t << 1));
    middle = ss_median3(Td, PA, middle - t, middle, middle + t);
    last   = ss_median3(Td, PA, last - 1 - (t << 1), last - 1 - t, last - 1);
    return ss_median3(Td, PA, first, middle, last);
}


/*---------------------------------------------------------------------------*/

/* Binary partition for substrings. */
static INLINE
int *
ss_partition(const int *PA,
                                        int *first, int *last, int depth) {
    int *a, *b;
    int t;
    for(a = first - 1, b = last;;) {
        for(; (++a < b) && ((PA[*a] + depth) >= (PA[*a + 1] + 1));) { *a = ~*a; }
        for(; (a < --b) && ((PA[*b] + depth) <  (PA[*b + 1] + 1));) { }
        if(b <= a) { break; }
        t = ~*b;
        *b = *a;
        *a = t;
    }
    if(first < a) { *first = ~*first; }
    return a;
}

/* Multikey introsort for medium size groups. */
static
void
ss_mintrosort(const unsigned char *T, const int *PA,
                            int *first, int *last,
                            int depth) {
#define STACK_SIZE SS_MISORT_STACKSIZE
    struct { int *a, *b, c; int d; } stack[STACK_SIZE];
    const unsigned char *Td;
    int *a, *b, *c, *d, *e, *f;
    int s, t;
    int ssize;
    int limit;
    int v, x = 0;

    for(ssize = 0, limit = ss_ilg(last - first);;) {

        if((last - first) <= SS_INSERTIONSORT_THRESHOLD) {
#if 1 < SS_INSERTIONSORT_THRESHOLD
            if(1 < (last - first)) { ss_insertionsort(T, PA, first, last, depth); }
#endif
            STACK_POP(first, last, depth, limit);
            continue;
        }

        Td = T + depth;
        if(limit-- == 0) { ss_heapsort(Td, PA, first, last - first); }
        if(limit < 0) {
            for(a = first + 1, v = Td[PA[*first]]; a < last; ++a) {
                if((x = Td[PA[*a]]) != v) {
                    if(1 < (a - first)) { break; }
                    v = x;
                    first = a;
                }
            }
            if(Td[PA[*first] - 1] < v) {
                first = ss_partition(PA, first, a, depth);
            }
            if((a - first) <= (last - a)) {
                if(1 < (a - first)) {
                    STACK_PUSH(a, last, depth, -1);
                    last = a, depth += 1, limit = ss_ilg(a - first);
                } else {
                    first = a, limit = -1;
                }
            } else {
                if(1 < (last - a)) {
                    STACK_PUSH(first, a, depth + 1, ss_ilg(a - first));
                    first = a, limit = -1;
                } else {
                    last = a, depth += 1, limit = ss_ilg(a - first);
                }
            }
            continue;
        }

        /* choose pivot */
        a = ss_pivot(Td, PA, first, last);
        v = Td[PA[*a]];
        SWAP(*first, *a);

        /* partition */
        for(b = first; (++b < last) && ((x = Td[PA[*b]]) == v);) { }
        if(((a = b) < last) && (x < v)) {
            for(; (++b < last) && ((x = Td[PA[*b]]) <= v);) {
                if(x == v) { SWAP(*b, *a); ++a; }
            }
        }
        for(c = last; (b < --c) && ((x = Td[PA[*c]]) == v);) { }
        if((b < (d = c)) && (x > v)) {
            for(; (b < --c) && ((x = Td[PA[*c]]) >= v);) {
                if(x == v) { SWAP(*c, *d); --d; }
            }
        }
        for(; b < c;) {
            SWAP(*b, *c);
            for(; (++b < c) && ((x = Td[PA[*b]]) <= v);) {
                if(x == v) { SWAP(*b, *a); ++a; }
            }
            for(; (b < --c) && ((x = Td[PA[*c]]) >= v);) {
                if(x == v) { SWAP(*c, *d); --d; }
            }
        }

        if(a <= d) {
            c = b - 1;

            if((s = a - first) > (t = b - a)) { s = t; }
            for(e = first, f = b - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }
            if((s = d - c) > (t = last - d - 1)) { s = t; }
            for(e = b, f = last - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }

            a = first + (b - a), c = last - (d - c);
            b = (v <= Td[PA[*a] - 1]) ? a : ss_partition(PA, a, c, depth);

            if((a - first) <= (last - c)) {
                if((last - c) <= (c - b)) {
                    STACK_PUSH(b, c, depth + 1, ss_ilg(c - b));
                    STACK_PUSH(c, last, depth, limit);
                    last = a;
                } else if((a - first) <= (c - b)) {
                    STACK_PUSH(c, last, depth, limit);
                    STACK_PUSH(b, c, depth + 1, ss_ilg(c - b));
                    last = a;
                } else {
                    STACK_PUSH(c, last, depth, limit);
                    STACK_PUSH(first, a, depth, limit);
                    first = b, last = c, depth += 1, limit = ss_ilg(c - b);
                }
            } else {
                if((a - first) <= (c - b)) {
                    STACK_PUSH(b, c, depth + 1, ss_ilg(c - b));
                    STACK_PUSH(first, a, depth, limit);
                    first = c;
                } else if((last - c) <= (c - b)) {
                    STACK_PUSH(first, a, depth, limit);
                    STACK_PUSH(b, c, depth + 1, ss_ilg(c - b));
                    first = c;
                } else {
                    STACK_PUSH(first, a, depth, limit);
                    STACK_PUSH(c, last, depth, limit);
                    first = b, last = c, depth += 1, limit = ss_ilg(c - b);
                }
            }
        } else {
            limit += 1;
            if(Td[PA[*first] - 1] < v) {
                first = ss_partition(PA, first, last, depth);
                limit = ss_ilg(last - first);
            }
            depth += 1;
        }
    }
#undef STACK_SIZE
}

#endif /* (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE) */


/*---------------------------------------------------------------------------*/

#if SS_BLOCKSIZE != 0

static INLINE
void
ss_blockswap(int *a, int *b, int n) {
    int t;
    for(; 0 < n; --n, ++a, ++b) {
        t = *a, *a = *b, *b = t;
    }
}

static INLINE
void
ss_rotate(int *first, int *middle, int *last) {
    int *a, *b, t;
    int l, r;
    l = middle - first, r = last - middle;
    for(; (0 < l) && (0 < r);) {
        if(l == r) { ss_blockswap(first, middle, l); break; }
        if(l < r) {
            a = last - 1, b = middle - 1;
            t = *a;
            do {
                *a-- = *b, *b-- = *a;
                if(b < first) {
                    *a = t;
                    last = a;
                    if((r -= l + 1) <= l) { break; }
                    a -= 1, b = middle - 1;
                    t = *a;
                }
            } while(1);
        } else {
            a = first, b = middle;
            t = *a;
            do {
                *a++ = *b, *b++ = *a;
                if(last <= b) {
                    *a = t;
                    first = a + 1;
                    if((l -= r + 1) <= r) { break; }
                    a += 1, b = middle;
                    t = *a;
                }
            } while(1);
        }
    }
}


/*---------------------------------------------------------------------------*/

static
void
ss_inplacemerge(const unsigned char *T, const int *PA,
                                int *first, int *middle, int *last,
                                int depth) {
    const int *p;
    int *a, *b;
    int len, half;
    int q, r;
    int x;

    for(;;) {
        if(*(last - 1) < 0) { x = 1; p = PA + ~*(last - 1); }
        else                { x = 0; p = PA +  *(last - 1); }
        for(a = first, len = middle - first, half = len >> 1, r = -1;
                0 < len;
                len = half, half >>= 1) {
            b = a + half;
            q = ss_compare(T, PA + ((0 <= *b) ? *b : ~*b), p, depth);
            if(q < 0) {
                a = b + 1;
                half -= (len & 1) ^ 1;
            } else {
                r = q;
            }
        }
        if(a < middle) {
            if(r == 0) { *a = ~*a; }
            ss_rotate(a, middle, last);
            last -= middle - a;
            middle = a;
            if(first == middle) { break; }
        }
        --last;
        if(x != 0) { while(*--last < 0) { } }
        if(middle == last) { break; }
    }
}


/*---------------------------------------------------------------------------*/

/* Merge-forward with internal buffer. */
static
void
ss_mergeforward(const unsigned char *T, const int *PA,
                                int *first, int *middle, int *last,
                                int *buf, int depth) {
    int *a, *b, *c, *bufend;
    int t;
    int r;

    bufend = buf + (middle - first) - 1;
    ss_blockswap(buf, first, middle - first);

    for(t = *(a = first), b = buf, c = middle;;) {
        r = ss_compare(T, PA + *b, PA + *c, depth);
        if(r < 0) {
            do {
                *a++ = *b;
                if(bufend <= b) { *bufend = t; return; }
                *b++ = *a;
            } while(*b < 0);
        } else if(r > 0) {
            do {
                *a++ = *c, *c++ = *a;
                if(last <= c) {
                    while(b < bufend) { *a++ = *b, *b++ = *a; }
                    *a = *b, *b = t;
                    return;
                }
            } while(*c < 0);
        } else {
            *c = ~*c;
            do {
                *a++ = *b;
                if(bufend <= b) { *bufend = t; return; }
                *b++ = *a;
            } while(*b < 0);

            do {
                *a++ = *c, *c++ = *a;
                if(last <= c) {
                    while(b < bufend) { *a++ = *b, *b++ = *a; }
                    *a = *b, *b = t;
                    return;
                }
            } while(*c < 0);
        }
    }
}

/* Merge-backward with internal buffer. */
static
void
ss_mergebackward(const unsigned char *T, const int *PA,
                                 int *first, int *middle, int *last,
                                 int *buf, int depth) {
    const int *p1, *p2;
    int *a, *b, *c, *bufend;
    int t;
    int r;
    int x;

    bufend = buf + (last - middle) - 1;
    ss_blockswap(buf, middle, last - middle);

    x = 0;
    if(*bufend < 0)       { p1 = PA + ~*bufend; x |= 1; }
    else                  { p1 = PA +  *bufend; }
    if(*(middle - 1) < 0) { p2 = PA + ~*(middle - 1); x |= 2; }
    else                  { p2 = PA +  *(middle - 1); }
    for(t = *(a = last - 1), b = bufend, c = middle - 1;;) {
        r = ss_compare(T, p1, p2, depth);
        if(0 < r) {
            if(x & 1) { do { *a-- = *b, *b-- = *a; } while(*b < 0); x ^= 1; }
            *a-- = *b;
            if(b <= buf) { *buf = t; break; }
            *b-- = *a;
            if(*b < 0) { p1 = PA + ~*b; x |= 1; }
            else       { p1 = PA +  *b; }
        } else if(r < 0) {
            if(x & 2) { do { *a-- = *c, *c-- = *a; } while(*c < 0); x ^= 2; }
            *a-- = *c, *c-- = *a;
            if(c < first) {
                while(buf < b) { *a-- = *b, *b-- = *a; }
                *a = *b, *b = t;
                break;
            }
            if(*c < 0) { p2 = PA + ~*c; x |= 2; }
            else       { p2 = PA +  *c; }
        } else {
            if(x & 1) { do { *a-- = *b, *b-- = *a; } while(*b < 0); x ^= 1; }
            *a-- = ~*b;
            if(b <= buf) { *buf = t; break; }
            *b-- = *a;
            if(x & 2) { do { *a-- = *c, *c-- = *a; } while(*c < 0); x ^= 2; }
            *a-- = *c, *c-- = *a;
            if(c < first) {
                while(buf < b) { *a-- = *b, *b-- = *a; }
                *a = *b, *b = t;
                break;
            }
            if(*b < 0) { p1 = PA + ~*b; x |= 1; }
            else       { p1 = PA +  *b; }
            if(*c < 0) { p2 = PA + ~*c; x |= 2; }
            else       { p2 = PA +  *c; }
        }
    }
}

/* D&C based merge. */
static
void
ss_swapmerge(const unsigned char *T, const int *PA,
                         int *first, int *middle, int *last,
                         int *buf, int bufsize, int depth) {
#define STACK_SIZE SS_SMERGE_STACKSIZE
#define GETIDX(a) ((0 <= (a)) ? (a) : (~(a)))
#define MERGE_CHECK(a, b, c)\
    do {\
        if(((c) & 1) ||\
             (((c) & 2) && (ss_compare(T, PA + GETIDX(*((a) - 1)), PA + *(a), depth) == 0))) {\
            *(a) = ~*(a);\
        }\
        if(((c) & 4) && ((ss_compare(T, PA + GETIDX(*((b) - 1)), PA + *(b), depth) == 0))) {\
            *(b) = ~*(b);\
        }\
    } while(0)
    struct { int *a, *b, *c; int d; } stack[STACK_SIZE];
    int *l, *r, *lm, *rm;
    int m, len, half;
    int ssize;
    int check, next;

    for(check = 0, ssize = 0;;) {
        if((last - middle) <= bufsize) {
            if((first < middle) && (middle < last)) {
                ss_mergebackward(T, PA, first, middle, last, buf, depth);
            }
            MERGE_CHECK(first, last, check);
            STACK_POP(first, middle, last, check);
            continue;
        }

        if((middle - first) <= bufsize) {
            if(first < middle) {
                ss_mergeforward(T, PA, first, middle, last, buf, depth);
            }
            MERGE_CHECK(first, last, check);
            STACK_POP(first, middle, last, check);
            continue;
        }

        for(m = 0, len = MIN(middle - first, last - middle), half = len >> 1;
                0 < len;
                len = half, half >>= 1) {
            if(ss_compare(T, PA + GETIDX(*(middle + m + half)),
                                             PA + GETIDX(*(middle - m - half - 1)), depth) < 0) {
                m += half + 1;
                half -= (len & 1) ^ 1;
            }
        }

        if(0 < m) {
            lm = middle - m, rm = middle + m;
            ss_blockswap(lm, middle, m);
            l = r = middle, next = 0;
            if(rm < last) {
                if(*rm < 0) {
                    *rm = ~*rm;
                    if(first < lm) { for(; *--l < 0;) { } next |= 4; }
                    next |= 1;
                } else if(first < lm) {
                    for(; *r < 0; ++r) { }
                    next |= 2;
                }
            }

            if((l - first) <= (last - r)) {
                STACK_PUSH(r, rm, last, (next & 3) | (check & 4));
                middle = lm, last = l, check = (check & 3) | (next & 4);
            } else {
                if((next & 2) && (r == middle)) { next ^= 6; }
                STACK_PUSH(first, lm, l, (check & 3) | (next & 4));
                first = r, middle = rm, check = (next & 3) | (check & 4);
            }
        } else {
            if(ss_compare(T, PA + GETIDX(*(middle - 1)), PA + *middle, depth) == 0) {
                *middle = ~*middle;
            }
            MERGE_CHECK(first, last, check);
            STACK_POP(first, middle, last, check);
        }
    }
#undef STACK_SIZE
}

#endif /* SS_BLOCKSIZE != 0 */


/*---------------------------------------------------------------------------*/

/* Substring sort */
static
void
sssort(const unsigned char *T, const int *PA,
             int *first, int *last,
             int *buf, int bufsize,
             int depth, int n, int lastsuffix) {
    int *a;
#if SS_BLOCKSIZE != 0
    int *b, *middle, *curbuf;
    int j, k, curbufsize, limit;
#endif
    int i;

    if(lastsuffix != 0) { ++first; }

#if SS_BLOCKSIZE == 0
    ss_mintrosort(T, PA, first, last, depth);
#else
    if((bufsize < SS_BLOCKSIZE) &&
            (bufsize < (last - first)) &&
            (bufsize < (limit = ss_isqrt(last - first)))) {
        if(SS_BLOCKSIZE < limit) { limit = SS_BLOCKSIZE; }
        buf = middle = last - limit, bufsize = limit;
    } else {
        middle = last, limit = 0;
    }
    for(a = first, i = 0; SS_BLOCKSIZE < (middle - a); a += SS_BLOCKSIZE, ++i) {
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
        ss_mintrosort(T, PA, a, a + SS_BLOCKSIZE, depth);
#elif 1 < SS_BLOCKSIZE
        ss_insertionsort(T, PA, a, a + SS_BLOCKSIZE, depth);
#endif
        curbufsize = last - (a + SS_BLOCKSIZE);
        curbuf = a + SS_BLOCKSIZE;
        if(curbufsize <= bufsize) { curbufsize = bufsize, curbuf = buf; }
        for(b = a, k = SS_BLOCKSIZE, j = i; j & 1; b -= k, k <<= 1, j >>= 1) {
            ss_swapmerge(T, PA, b - k, b, b + k, curbuf, curbufsize, depth);
        }
    }
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
    ss_mintrosort(T, PA, a, middle, depth);
#elif 1 < SS_BLOCKSIZE
    ss_insertionsort(T, PA, a, middle, depth);
#endif
    for(k = SS_BLOCKSIZE; i != 0; k <<= 1, i >>= 1) {
        if(i & 1) {
            ss_swapmerge(T, PA, a - k, a, middle, buf, bufsize, depth);
            a -= k;
        }
    }
    if(limit != 0) {
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
        ss_mintrosort(T, PA, middle, last, depth);
#elif 1 < SS_BLOCKSIZE
        ss_insertionsort(T, PA, middle, last, depth);
#endif
        ss_inplacemerge(T, PA, first, middle, last, depth);
    }
#endif

    if(lastsuffix != 0) {
        /* Insert last type B* suffix. */
        int PAi[2]; PAi[0] = PA[*(first - 1)], PAi[1] = n - 2;
        for(a = first, i = *(first - 1);
                (a < last) && ((*a < 0) || (0 < ss_compare(T, &(PAi[0]), PA + *a, depth)));
                ++a) {
            *(a - 1) = *a;
        }
        *(a - 1) = i;
    }
}


/*---------------------------------------------------------------------------*/

static INLINE
int
tr_ilg(int n) {
    return (n & 0xffff0000) ?
                    ((n & 0xff000000) ?
                        24 + lg_table[(n >> 24) & 0xff] :
                        16 + lg_table[(n >> 16) & 0xff]) :
                    ((n & 0x0000ff00) ?
                         8 + lg_table[(n >>  8) & 0xff] :
                         0 + lg_table[(n >>  0) & 0xff]);
}


/*---------------------------------------------------------------------------*/

/* Simple insertionsort for small size groups. */
static
void
tr_insertionsort(const int *ISAd, int *first, int *last) {
    int *a, *b;
    int t, r;

    for(a = first + 1; a < last; ++a) {
        for(t = *a, b = a - 1; 0 > (r = ISAd[t] - ISAd[*b]);) {
            do { *(b + 1) = *b; } while((first <= --b) && (*b < 0));
            if(b < first) { break; }
        }
        if(r == 0) { *b = ~*b; }
        *(b + 1) = t;
    }
}


/*---------------------------------------------------------------------------*/

static INLINE
void
tr_fixdown(const int *ISAd, int *SA, int i, int size) {
    int j, k;
    int v;
    int c, d, e;

    for(v = SA[i], c = ISAd[v]; (j = 2 * i + 1) < size; SA[i] = SA[k], i = k) {
        d = ISAd[SA[k = j++]];
        if(d < (e = ISAd[SA[j]])) { k = j; d = e; }
        if(d <= c) { break; }
    }
    SA[i] = v;
}

/* Simple top-down heapsort. */
static
void
tr_heapsort(const int *ISAd, int *SA, int size) {
    int i, m;
    int t;

    m = size;
    if((size % 2) == 0) {
        m--;
        if(ISAd[SA[m / 2]] < ISAd[SA[m]]) { SWAP(SA[m], SA[m / 2]); }
    }

    for(i = m / 2 - 1; 0 <= i; --i) { tr_fixdown(ISAd, SA, i, m); }
    if((size % 2) == 0) { SWAP(SA[0], SA[m]); tr_fixdown(ISAd, SA, 0, m); }
    for(i = m - 1; 0 < i; --i) {
        t = SA[0], SA[0] = SA[i];
        tr_fixdown(ISAd, SA, 0, i);
        SA[i] = t;
    }
}


/*---------------------------------------------------------------------------*/

/* Returns the median of three elements. */
static INLINE
int *
tr_median3(const int *ISAd, int *v1, int *v2, int *v3) {
    int *t;
    if(ISAd[*v1] > ISAd[*v2]) { SWAP(v1, v2); }
    if(ISAd[*v2] > ISAd[*v3]) {
        if(ISAd[*v1] > ISAd[*v3]) { return v1; }
        else { return v3; }
    }
    return v2;
}

/* Returns the median of five elements. */
static INLINE
int *
tr_median5(const int *ISAd,
                     int *v1, int *v2, int *v3, int *v4, int *v5) {
    int *t;
    if(ISAd[*v2] > ISAd[*v3]) { SWAP(v2, v3); }
    if(ISAd[*v4] > ISAd[*v5]) { SWAP(v4, v5); }
    if(ISAd[*v2] > ISAd[*v4]) { SWAP(v2, v4); SWAP(v3, v5); }
    if(ISAd[*v1] > ISAd[*v3]) { SWAP(v1, v3); }
    if(ISAd[*v1] > ISAd[*v4]) { SWAP(v1, v4); SWAP(v3, v5); }
    if(ISAd[*v3] > ISAd[*v4]) { return v4; }
    return v3;
}

/* Returns the pivot element. */
static INLINE
int *
tr_pivot(const int *ISAd, int *first, int *last) {
    int *middle;
    int t;

    t = last - first;
    middle = first + t / 2;

    if(t <= 512) {
        if(t <= 32) {
            return tr_median3(ISAd, first, middle, last - 1);
        } else {
            t >>= 2;
            return tr_median5(ISAd, first, first + t, middle, last - 1 - t, last - 1);
        }
    }
    t >>= 3;
    first  = tr_median3(ISAd, first, first + t, first + (t << 1));
    middle = tr_median3(ISAd, middle - t, middle, middle + t);
    last   = tr_median3(ISAd, last - 1 - (t << 1), last - 1 - t, last - 1);
    return tr_median3(ISAd, first, middle, last);
}


/*---------------------------------------------------------------------------*/

typedef struct _trbudget_t trbudget_t;
struct _trbudget_t {
    int chance;
    int remain;
    int incval;
    int count;
};

static INLINE
void
trbudget_init(trbudget_t *budget, int chance, int incval) {
    budget->chance = chance;
    budget->remain = budget->incval = incval;
}

static INLINE
int
trbudget_check(trbudget_t *budget, int size) {
    if(size <= budget->remain) { budget->remain -= size; return 1; }
    if(budget->chance == 0) { budget->count += size; return 0; }
    budget->remain += budget->incval - size;
    budget->chance -= 1;
    return 1;
}


/*---------------------------------------------------------------------------*/

static INLINE
void
tr_partition(const int *ISAd,
                         int *first, int *middle, int *last,
                         int **pa, int **pb, int v) {
    int *a, *b, *c, *d, *e, *f;
    int t, s;
    int x = 0;

    for(b = middle - 1; (++b < last) && ((x = ISAd[*b]) == v);) { }
    if(((a = b) < last) && (x < v)) {
        for(; (++b < last) && ((x = ISAd[*b]) <= v);) {
            if(x == v) { SWAP(*b, *a); ++a; }
        }
    }
    for(c = last; (b < --c) && ((x = ISAd[*c]) == v);) { }
    if((b < (d = c)) && (x > v)) {
        for(; (b < --c) && ((x = ISAd[*c]) >= v);) {
            if(x == v) { SWAP(*c, *d); --d; }
        }
    }
    for(; b < c;) {
        SWAP(*b, *c);
        for(; (++b < c) && ((x = ISAd[*b]) <= v);) {
            if(x == v) { SWAP(*b, *a); ++a; }
        }
        for(; (b < --c) && ((x = ISAd[*c]) >= v);) {
            if(x == v) { SWAP(*c, *d); --d; }
        }
    }

    if(a <= d) {
        c = b - 1;
        if((s = a - first) > (t = b - a)) { s = t; }
        for(e = first, f = b - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }
        if((s = d - c) > (t = last - d - 1)) { s = t; }
        for(e = b, f = last - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }
        first += (b - a), last -= (d - c);
    }
    *pa = first, *pb = last;
}

static
void
tr_copy(int *ISA, const int *SA,
                int *first, int *a, int *b, int *last,
                int depth) {
    /* sort suffixes of middle partition
         by using sorted order of suffixes of left and right partition. */
    int *c, *d, *e;
    int s, v;

    v = b - SA - 1;
    for(c = first, d = a - 1; c <= d; ++c) {
        if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
            *++d = s;
            ISA[s] = d - SA;
        }
    }
    for(c = last - 1, e = d + 1, d = b; e < d; --c) {
        if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
            *--d = s;
            ISA[s] = d - SA;
        }
    }
}

static
void
tr_partialcopy(int *ISA, const int *SA,
                             int *first, int *a, int *b, int *last,
                             int depth) {
    int *c, *d, *e;
    int s, v;
    int rank, lastrank, newrank = -1;

    v = b - SA - 1;
    lastrank = -1;
    for(c = first, d = a - 1; c <= d; ++c) {
        if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
            *++d = s;
            rank = ISA[s + depth];
            if(lastrank != rank) { lastrank = rank; newrank = d - SA; }
            ISA[s] = newrank;
        }
    }

    lastrank = -1;
    for(e = d; first <= e; --e) {
        rank = ISA[*e];
        if(lastrank != rank) { lastrank = rank; newrank = e - SA; }
        if(newrank != rank) { ISA[*e] = newrank; }
    }

    lastrank = -1;
    for(c = last - 1, e = d + 1, d = b; e < d; --c) {
        if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
            *--d = s;
            rank = ISA[s + depth];
            if(lastrank != rank) { lastrank = rank; newrank = d - SA; }
            ISA[s] = newrank;
        }
    }
}

static
void
tr_introsort(int *ISA, const int *ISAd,
                         int *SA, int *first, int *last,
                         trbudget_t *budget) {
#define STACK_SIZE TR_STACKSIZE
    struct { const int *a; int *b, *c; int d, e; }stack[STACK_SIZE];
    int *a, *b, *c;
    int t;
    int v, x = 0;
    int incr = ISAd - ISA;
    int limit, next;
    int ssize, trlink = -1;

    for(ssize = 0, limit = tr_ilg(last - first);;) {

        if(limit < 0) {
            if(limit == -1) {
                /* tandem repeat partition */
                tr_partition(ISAd - incr, first, first, last, &a, &b, last - SA - 1);

                /* update ranks */
                if(a < last) {
                    for(c = first, v = a - SA - 1; c < a; ++c) { ISA[*c] = v; }
                }
                if(b < last) {
                    for(c = a, v = b - SA - 1; c < b; ++c) { ISA[*c] = v; }
                }

                /* push */
                if(1 < (b - a)) {
                    STACK_PUSH5(NULL, a, b, 0, 0);
                    STACK_PUSH5(ISAd - incr, first, last, -2, trlink);
                    trlink = ssize - 2;
                }
                if((a - first) <= (last - b)) {
                    if(1 < (a - first)) {
                        STACK_PUSH5(ISAd, b, last, tr_ilg(last - b), trlink);
                        last = a, limit = tr_ilg(a - first);
                    } else if(1 < (last - b)) {
                        first = b, limit = tr_ilg(last - b);
                    } else {
                        STACK_POP5(ISAd, first, last, limit, trlink);
                    }
                } else {
                    if(1 < (last - b)) {
                        STACK_PUSH5(ISAd, first, a, tr_ilg(a - first), trlink);
                        first = b, limit = tr_ilg(last - b);
                    } else if(1 < (a - first)) {
                        last = a, limit = tr_ilg(a - first);
                    } else {
                        STACK_POP5(ISAd, first, last, limit, trlink);
                    }
                }
            } else if(limit == -2) {
                /* tandem repeat copy */
                a = stack[--ssize].b, b = stack[ssize].c;
                if(stack[ssize].d == 0) {
                    tr_copy(ISA, SA, first, a, b, last, ISAd - ISA);
                } else {
                    if(0 <= trlink) { stack[trlink].d = -1; }
                    tr_partialcopy(ISA, SA, first, a, b, last, ISAd - ISA);
                }
                STACK_POP5(ISAd, first, last, limit, trlink);
            } else {
                /* sorted partition */
                if(0 <= *first) {
                    a = first;
                    do { ISA[*a] = a - SA; } while((++a < last) && (0 <= *a));
                    first = a;
                }
                if(first < last) {
                    a = first; do { *a = ~*a; } while(*++a < 0);
                    next = (ISA[*a] != ISAd[*a]) ? tr_ilg(a - first + 1) : -1;
                    if(++a < last) { for(b = first, v = a - SA - 1; b < a; ++b) { ISA[*b] = v; } }

                    /* push */
                    if(trbudget_check(budget, a - first)) {
                        if((a - first) <= (last - a)) {
                            STACK_PUSH5(ISAd, a, last, -3, trlink);
                            ISAd += incr, last = a, limit = next;
                        } else {
                            if(1 < (last - a)) {
                                STACK_PUSH5(ISAd + incr, first, a, next, trlink);
                                first = a, limit = -3;
                            } else {
                                ISAd += incr, last = a, limit = next;
                            }
                        }
                    } else {
                        if(0 <= trlink) { stack[trlink].d = -1; }
                        if(1 < (last - a)) {
                            first = a, limit = -3;
                        } else {
                            STACK_POP5(ISAd, first, last, limit, trlink);
                        }
                    }
                } else {
                    STACK_POP5(ISAd, first, last, limit, trlink);
                }
            }
            continue;
        }

        if((last - first) <= TR_INSERTIONSORT_THRESHOLD) {
            tr_insertionsort(ISAd, first, last);
            limit = -3;
            continue;
        }

        if(limit-- == 0) {
            tr_heapsort(ISAd, first, last - first);
            for(a = last - 1; first < a; a = b) {
                for(x = ISAd[*a], b = a - 1; (first <= b) && (ISAd[*b] == x); --b) { *b = ~*b; }
            }
            limit = -3;
            continue;
        }

        /* choose pivot */
        a = tr_pivot(ISAd, first, last);
        SWAP(*first, *a);
        v = ISAd[*first];

        /* partition */
        tr_partition(ISAd, first, first + 1, last, &a, &b, v);
        if((last - first) != (b - a)) {
            next = (ISA[*a] != v) ? tr_ilg(b - a) : -1;

            /* update ranks */
            for(c = first, v = a - SA - 1; c < a; ++c) { ISA[*c] = v; }
            if(b < last) { for(c = a, v = b - SA - 1; c < b; ++c) { ISA[*c] = v; } }

            /* push */
            if((1 < (b - a)) && (trbudget_check(budget, b - a))) {
                if((a - first) <= (last - b)) {
                    if((last - b) <= (b - a)) {
                        if(1 < (a - first)) {
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            STACK_PUSH5(ISAd, b, last, limit, trlink);
                            last = a;
                        } else if(1 < (last - b)) {
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            first = b;
                        } else {
                            ISAd += incr, first = a, last = b, limit = next;
                        }
                    } else if((a - first) <= (b - a)) {
                        if(1 < (a - first)) {
                            STACK_PUSH5(ISAd, b, last, limit, trlink);
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            last = a;
                        } else {
                            STACK_PUSH5(ISAd, b, last, limit, trlink);
                            ISAd += incr, first = a, last = b, limit = next;
                        }
                    } else {
                        STACK_PUSH5(ISAd, b, last, limit, trlink);
                        STACK_PUSH5(ISAd, first, a, limit, trlink);
                        ISAd += incr, first = a, last = b, limit = next;
                    }
                } else {
                    if((a - first) <= (b - a)) {
                        if(1 < (last - b)) {
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            STACK_PUSH5(ISAd, first, a, limit, trlink);
                            first = b;
                        } else if(1 < (a - first)) {
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            last = a;
                        } else {
                            ISAd += incr, first = a, last = b, limit = next;
                        }
                    } else if((last - b) <= (b - a)) {
                        if(1 < (last - b)) {
                            STACK_PUSH5(ISAd, first, a, limit, trlink);
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            first = b;
                        } else {
                            STACK_PUSH5(ISAd, first, a, limit, trlink);
                            ISAd += incr, first = a, last = b, limit = next;
                        }
                    } else {
                        STACK_PUSH5(ISAd, first, a, limit, trlink);
                        STACK_PUSH5(ISAd, b, last, limit, trlink);
                        ISAd += incr, first = a, last = b, limit = next;
                    }
                }
            } else {
                if((1 < (b - a)) && (0 <= trlink)) { stack[trlink].d = -1; }
                if((a - first) <= (last - b)) {
                    if(1 < (a - first)) {
                        STACK_PUSH5(ISAd, b, last, limit, trlink);
                        last = a;
                    } else if(1 < (last - b)) {
                        first = b;
                    } else {
                        STACK_POP5(ISAd, first, last, limit, trlink);
                    }
                } else {
                    if(1 < (last - b)) {
                        STACK_PUSH5(ISAd, first, a, limit, trlink);
                        first = b;
                    } else if(1 < (a - first)) {
                        last = a;
                    } else {
                        STACK_POP5(ISAd, first, last, limit, trlink);
                    }
                }
            }
        } else {
            if(trbudget_check(budget, last - first)) {
                limit = tr_ilg(last - first), ISAd += incr;
            } else {
                if(0 <= trlink) { stack[trlink].d = -1; }
                STACK_POP5(ISAd, first, last, limit, trlink);
            }
        }
    }
#undef STACK_SIZE
}



/*---------------------------------------------------------------------------*/

/* Tandem repeat sort */
static
void
trsort(int *ISA, int *SA, int n, int depth) {
    int *ISAd;
    int *first, *last;
    trbudget_t budget;
    int t, skip, unsorted;

    trbudget_init(&budget, tr_ilg(n) * 2 / 3, n);
/*  trbudget_init(&budget, tr_ilg(n) * 3 / 4, n); */
    for(ISAd = ISA + depth; -n < *SA; ISAd += ISAd - ISA) {
        first = SA;
        skip = 0;
        unsorted = 0;
        do {
            if((t = *first) < 0) { first -= t; skip += t; }
            else {
                if(skip != 0) { *(first + skip) = skip; skip = 0; }
                last = SA + ISA[t] + 1;
                if(1 < (last - first)) {
                    budget.count = 0;
                    tr_introsort(ISA, ISAd, SA, first, last, &budget);
                    if(budget.count != 0) { unsorted += budget.count; }
                    else { skip = first - last; }
                } else if((last - first) == 1) {
                    skip = -1;
                }
                first = last;
            }
        } while(first < (SA + n));
        if(skip != 0) { *(first + skip) = skip; }
        if(unsorted == 0) { break; }
    }
}


/*---------------------------------------------------------------------------*/

/* Sorts suffixes of type B*. */
static
int
sort_typeBstar(const unsigned char *T, int *SA,
                             int *bucket_A, int *bucket_B,
                             int n) {
    int *PAb, *ISAb, *buf;
#ifdef _OPENMP
    int *curbuf;
    int l;
#endif
    int i, j, k, t, m, bufsize;
    int c0, c1;
#ifdef _OPENMP
    int d0, d1;
    int tmp;
#endif

    /* Initialize bucket arrays. */
    for(i = 0; i < BUCKET_A_SIZE; ++i) { bucket_A[i] = 0; }
    for(i = 0; i < BUCKET_B_SIZE; ++i) { bucket_B[i] = 0; }

    /* Count the number of occurrences of the first one or two characters of each
         type A, B and B* suffix. Moreover, store the beginning position of all
         type B* suffixes into the array SA. */
    for(i = n - 1, m = n, c0 = T[n - 1]; 0 <= i;) {
        /* type A suffix. */
        do { ++BUCKET_A(c1 = c0); } while((0 <= --i) && ((c0 = T[i]) >= c1));
        if(0 <= i) {
            /* type B* suffix. */
            ++BUCKET_BSTAR(c0, c1);
            SA[--m] = i;
            /* type B suffix. */
            for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) {
                ++BUCKET_B(c0, c1);
            }
        }
    }
    m = n - m;
/*
note:
    A type B* suffix is lexicographically smaller than a type B suffix that
    begins with the same first two characters.
*/

    /* Calculate the index of start/end point of each bucket. */
    for(c0 = 0, i = 0, j = 0; c0 < ALPHABET_SIZE; ++c0) {
        t = i + BUCKET_A(c0);
        BUCKET_A(c0) = i + j; /* start point */
        i = t + BUCKET_B(c0, c0);
        for(c1 = c0 + 1; c1 < ALPHABET_SIZE; ++c1) {
            j += BUCKET_BSTAR(c0, c1);
            BUCKET_BSTAR(c0, c1) = j; /* end point */
            i += BUCKET_B(c0, c1);
        }
    }

    if(0 < m) {
        /* Sort the type B* suffixes by their first two characters. */
        PAb = SA + n - m; ISAb = SA + m;
        for(i = m - 2; 0 <= i; --i) {
            t = PAb[i], c0 = T[t], c1 = T[t + 1];
            SA[--BUCKET_BSTAR(c0, c1)] = i;
        }
        t = PAb[m - 1], c0 = T[t], c1 = T[t + 1];
        SA[--BUCKET_BSTAR(c0, c1)] = m - 1;

        /* Sort the type B* substrings using sssort. */
#ifdef _OPENMP
        tmp = omp_get_max_threads();
        buf = SA + m, bufsize = (n - (2 * m)) / tmp;
        c0 = ALPHABET_SIZE - 2, c1 = ALPHABET_SIZE - 1, j = m;
#pragma omp parallel default(shared) private(curbuf, k, l, d0, d1, tmp)
        {
            tmp = omp_get_thread_num();
            curbuf = buf + tmp * bufsize;
            k = 0;
            for(;;) {
                #pragma omp critical(sssort_lock)
                {
                    if(0 < (l = j)) {
                        d0 = c0, d1 = c1;
                        do {
                            k = BUCKET_BSTAR(d0, d1);
                            if(--d1 <= d0) {
                                d1 = ALPHABET_SIZE - 1;
                                if(--d0 < 0) { break; }
                            }
                        } while(((l - k) <= 1) && (0 < (l = k)));
                        c0 = d0, c1 = d1, j = k;
                    }
                }
                if(l == 0) { break; }
                sssort(T, PAb, SA + k, SA + l,
                             curbuf, bufsize, 2, n, *(SA + k) == (m - 1));
            }
        }
#else
        buf = SA + m, bufsize = n - (2 * m);
        for(c0 = ALPHABET_SIZE - 2, j = m; 0 < j; --c0) {
            for(c1 = ALPHABET_SIZE - 1; c0 < c1; j = i, --c1) {
                i = BUCKET_BSTAR(c0, c1);
                if(1 < (j - i)) {
                    sssort(T, PAb, SA + i, SA + j,
                                 buf, bufsize, 2, n, *(SA + i) == (m - 1));
                }
            }
        }
#endif

        /* Compute ranks of type B* substrings. */
        for(i = m - 1; 0 <= i; --i) {
            if(0 <= SA[i]) {
                j = i;
                do { ISAb[SA[i]] = i; } while((0 <= --i) && (0 <= SA[i]));
                SA[i + 1] = i - j;
                if(i <= 0) { break; }
            }
            j = i;
            do { ISAb[SA[i] = ~SA[i]] = j; } while(SA[--i] < 0);
            ISAb[SA[i]] = j;
        }

        /* Construct the inverse suffix array of type B* suffixes using trsort. */
        trsort(ISAb, SA, m, 1);

        /* Set the sorted order of tyoe B* suffixes. */
        for(i = n - 1, j = m, c0 = T[n - 1]; 0 <= i;) {
            for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) >= c1); --i, c1 = c0) { }
            if(0 <= i) {
                t = i;
                for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) { }
                SA[ISAb[--j]] = ((t == 0) || (1 < (t - i))) ? t : ~t;
            }
        }

        /* Calculate the index of start/end point of each bucket. */
        BUCKET_B(ALPHABET_SIZE - 1, ALPHABET_SIZE - 1) = n; /* end point */
        for(c0 = ALPHABET_SIZE - 2, k = m - 1; 0 <= c0; --c0) {
            i = BUCKET_A(c0 + 1) - 1;
            for(c1 = ALPHABET_SIZE - 1; c0 < c1; --c1) {
                t = i - BUCKET_B(c0, c1);
                BUCKET_B(c0, c1) = i; /* end point */

                /* Move all type B* suffixes to the correct position. */
                for(i = t, j = BUCKET_BSTAR(c0, c1);
                        j <= k;
                        --i, --k) { SA[i] = SA[k]; }
            }
            BUCKET_BSTAR(c0, c0 + 1) = i - BUCKET_B(c0, c0) + 1; /* start point */
            BUCKET_B(c0, c0) = i; /* end point */
        }
    }

    return m;
}

/* Constructs the suffix array by using the sorted order of type B* suffixes. */
static
void
construct_SA(const unsigned char *T, int *SA,
                         int *bucket_A, int *bucket_B,
                         int n, int m) {
    int *i, *j, *k;
    int s;
    int c0, c1, c2;

    if(0 < m) {
        /* Construct the sorted order of type B suffixes by using
             the sorted order of type B* suffixes. */
        for(c1 = ALPHABET_SIZE - 2; 0 <= c1; --c1) {
            /* Scan the suffix array from right to left. */
            for(i = SA + BUCKET_BSTAR(c1, c1 + 1),
                    j = SA + BUCKET_A(c1 + 1) - 1, k = NULL, c2 = -1;
                    i <= j;
                    --j) {
                if(0 < (s = *j)) {
                    assert(T[s] == c1);
                    assert(((s + 1) < n) && (T[s] <= T[s + 1]));
                    assert(T[s - 1] <= T[s]);
                    *j = ~s;
                    c0 = T[--s];
                    if((0 < s) && (T[s - 1] > c0)) { s = ~s; }
                    if(c0 != c2) {
                        if(0 <= c2) { BUCKET_B(c2, c1) = k - SA; }
                        k = SA + BUCKET_B(c2 = c0, c1);
                    }
                    assert(k < j);
                    *k-- = s;
                } else {
                    assert(((s == 0) && (T[s] == c1)) || (s < 0));
                    *j = ~s;
                }
            }
        }
    }

    /* Construct the suffix array by using
         the sorted order of type B suffixes. */
    k = SA + BUCKET_A(c2 = T[n - 1]);
    *k++ = (T[n - 2] < c2) ? ~(n - 1) : (n - 1);
    /* Scan the suffix array from left to right. */
    for(i = SA, j = SA + n; i < j; ++i) {
        if(0 < (s = *i)) {
            assert(T[s - 1] >= T[s]);
            c0 = T[--s];
            if((s == 0) || (T[s - 1] < c0)) { s = ~s; }
            if(c0 != c2) {
                BUCKET_A(c2) = k - SA;
                k = SA + BUCKET_A(c2 = c0);
            }
            assert(i < k);
            *k++ = s;
        } else {
            assert(s < 0);
            *i = ~s;
        }
    }
}

/* Constructs the burrows-wheeler transformed string directly
     by using the sorted order of type B* suffixes. */
static
int
construct_BWT(const unsigned char *T, int *SA,
                            int *bucket_A, int *bucket_B,
                            int n, int m) {
    int *i, *j, *k, *orig;
    int s;
    int c0, c1, c2;

    if(0 < m) {
        /* Construct the sorted order of type B suffixes by using
             the sorted order of type B* suffixes. */
        for(c1 = ALPHABET_SIZE - 2; 0 <= c1; --c1) {
            /* Scan the suffix array from right to left. */
            for(i = SA + BUCKET_BSTAR(c1, c1 + 1),
                    j = SA + BUCKET_A(c1 + 1) - 1, k = NULL, c2 = -1;
                    i <= j;
                    --j) {
                if(0 < (s = *j)) {
                    assert(T[s] == c1);
                    assert(((s + 1) < n) && (T[s] <= T[s + 1]));
                    assert(T[s - 1] <= T[s]);
                    c0 = T[--s];
                    *j = ~((int)c0);
                    if((0 < s) && (T[s - 1] > c0)) { s = ~s; }
                    if(c0 != c2) {
                        if(0 <= c2) { BUCKET_B(c2, c1) = k - SA; }
                        k = SA + BUCKET_B(c2 = c0, c1);
                    }
                    assert(k < j);
                    *k-- = s;
                } else if(s != 0) {
                    *j = ~s;
#ifndef NDEBUG
                } else {
                    assert(T[s] == c1);
#endif
                }
            }
        }
    }

    /* Construct the BWTed string by using
         the sorted order of type B suffixes. */
    k = SA + BUCKET_A(c2 = T[n - 1]);
    *k++ = (T[n - 2] < c2) ? ~((int)T[n - 2]) : (n - 1);
    /* Scan the suffix array from left to right. */
    for(i = SA, j = SA + n, orig = SA; i < j; ++i) {
        if(0 < (s = *i)) {
            assert(T[s - 1] >= T[s]);
            c0 = T[--s];
            *i = c0;
            if((0 < s) && (T[s - 1] < c0)) { s = ~((int)T[s - 1]); }
            if(c0 != c2) {
                BUCKET_A(c2) = k - SA;
                k = SA + BUCKET_A(c2 = c0);
            }
            assert(i < k);
            *k++ = s;
        } else if(s != 0) {
            *i = ~s;
        } else {
            orig = i;
        }
    }

    return orig - SA;
}


/*---------------------------------------------------------------------------*/

/*- Function -*/

int
bcm_divsufsort(const unsigned char *T, int *SA, int n) {
    int *bucket_A, *bucket_B;
    int m;
    int err = 0;

    /* Check arguments. */
    if((T == NULL) || (SA == NULL) || (n < 0)) { return -1; }
    else if(n == 0) { return 0; }
    else if(n == 1) { SA[0] = 0; return 0; }
    else if(n == 2) { m = (T[0] < T[1]); SA[m ^ 1] = 0, SA[m] = 1; return 0; }

    bucket_A = (int *)malloc(BUCKET_A_SIZE * sizeof(int));
    bucket_B = (int *)malloc(BUCKET_B_SIZE * sizeof(int));

    /* Suffixsort. */
    if((bucket_A != NULL) && (bucket_B != NULL)) {
        m = sort_typeBstar(T, SA, bucket_A, bucket_B, n);
        construct_SA(T, SA, bucket_A, bucket_B, n, m);
    } else {
        err = -2;
    }

    free(bucket_B);
    free(bucket_A);

    return err;
}

int
bcm_divbwt(const unsigned char *T, unsigned char *U, int *A, int n) {
    int *B;
    int *bucket_A, *bucket_B;
    int m, pidx, i;

    /* Check arguments. */
    if((T == NULL) || (U == NULL) || (n < 0)) { return -1; }
    else if(n <= 1) { if(n == 1) { U[0] = T[0]; } return n; }

    if((B = A) == NULL) { B = (int *)malloc((size_t)(n + 1) * sizeof(int)); }
    bucket_A = (int *)malloc(BUCKET_A_SIZE * sizeof(int));
    bucket_B = (int *)malloc(BUCKET_B_SIZE * sizeof(int));

    /* Burrows-Wheeler Transform. */
    if((B != NULL) && (bucket_A != NULL) && (bucket_B != NULL)) {
        m = sort_typeBstar(T, B, bucket_A, bucket_B, n);
        pidx = construct_BWT(T, B, bucket_A, bucket_B, n, m);

        /* Copy to output string. */
        U[0] = T[n - 1];
        for(i = 0; i < pidx; ++i) { U[i + 1] = (unsigned char)B[i]; }
        for(i += 1; i < n; ++i) { U[i] = (unsigned char)B[i]; }
        pidx += 1;
    } else {
        pidx = -2;
    }

    free(bucket_B);
    free(bucket_A);
    if(A == NULL) { free(B); }

    return pidx;
}

#endif // BCM_C

#line 1 "amalgamated_bcm.c" 
// BCM 1.40 - A BWT-based file compressor
// Written and placed in the public domain by Ilya Muravyov (UNLICENSE)
// Additional code by @r-lyeh (UNLICENSE)
//
// Notes:
// - BCM decoder has no dependencies.
// - BCM encoder requires libdivsufsort, which is MIT licensed.
// - #define BCM_NO_ENCODER if you want to exclude libdivsufsort from linkage.

unsigned bcm_bounds(unsigned inlen, unsigned flags);
unsigned bcm_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags/*[0..(4)..9]*/);
unsigned bcm_decode(const void *in, unsigned inlen, void *out, unsigned outlen);

// ---

#ifdef BCM_C
#pragma once
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if INTPTR_MAX >= INT64_MAX
#define BCM_64BITS 1
#else
#define BCM_64BITS 0
#endif

#ifndef BCM_REALLOC
#define BCM_REALLOC realloc
#endif

#  if defined _MSC_VER && !defined __thread
#define __thread __declspec(thread)
#elif defined __TINYC__ && !defined __thread
#define __thread __declspec(thread)
#endif

#ifndef BALZ_C
typedef struct mfile {
        uint8_t *begin, *seek, *end;
} mfile;
int minit(mfile *f, const void *ptr, int len) {
        f->begin = f->seek = f->end = (uint8_t*)ptr;
        f->end += len;
        return 0;
}
int mread(mfile *m, void *buf, int len) {
        if( len >= (m->end - m->seek) ) len = (m->end - m->seek);
        memcpy(buf,m->seek,len); m->seek += len;
        return len;
}
int mwrite(mfile *m, const void *buf, int len) {
        if( len >= (m->end - m->seek) ) len = (m->end - m->seek);
        memcpy(m->seek,buf,len); m->seek += len;
        return len;
}
int mtell(mfile *m) {
        return m->seek - m->begin;
}
int mavail(mfile *m) {
        return m->end - m->seek;
}
int mputc(mfile *m, int i) {
        uint8_t ch = i;
        return mwrite(m, &ch, 1);
}
int mgetc(mfile *m) { 
        if( mavail(m) <= 0 ) return -1;
        uint8_t ch; mread(m, &ch, 1); return ch;
}
#endif

int bcm_divbwt(const unsigned char *T, unsigned char *U, int *A, int n);

// Globals

static __thread mfile* g_in;
static __thread mfile* g_out;

typedef struct bcmEncode
{
    uint32_t low;
    uint32_t high;
    uint32_t code;
} bcmEncoder;

    void bcmeCtor(bcmEncoder *e)
    {
        e->low=0;
        e->high=0xFFFFFFFF;
        e->code=0;
    }

    void bcmeFlush(bcmEncoder *e)
    {
        for (int i=0; i<4; ++i)
        {
            mputc(g_out, e->low>>24);
            e->low<<=8;
        }
    }

    void bcmeInit(bcmEncoder *e)
    {
        for (int i=0; i<4; ++i)
            e->code=(e->code<<8)+mgetc(g_in);
    }

    void bcmeEncodeDirectBits(bcmEncoder *e, int N, uint32_t x)
    {
        for (uint32_t i=1<<(N-1); i!=0; i>>=1)
        {
            if (x&i)
                e->high=e->low+((e->high-e->low)>>1);
            else
                e->low+=((e->high-e->low)>>1)+1;

            if ((e->low^e->high)<(1<<24))
            {
                mputc(g_out, e->low>>24);
                e->low<<=8;
                e->high=(e->high<<8)+255;
            }
        }
    }

    void bcmeEncodeBit1(bcmEncoder *e, uint32_t p)
    {
#if BCM_64BITS
        e->high=e->low+(((uint64_t)(e->high-e->low)*p)>>18);
#else
        e->high=e->low+(((uint64_t)(e->high-e->low)*(p<<(32-18)))>>32);
#endif
        while ((e->low^e->high)<(1<<24))
        {
            mputc(g_out, e->low>>24);
            e->low<<=8;
            e->high=(e->high<<8)+255;
        }
    }

    void bcmeEncodeBit0(bcmEncoder *e, uint32_t p)
    {
#if BCM_64BITS
        e->low+=(((uint64_t)(e->high-e->low)*p)>>18)+1;
#else
        e->low+=(((uint64_t)(e->high-e->low)*(p<<(32-18)))>>32)+1;
#endif
        while ((e->low^e->high)<(1<<24))
        {
            mputc(g_out, e->low>>24);
            e->low<<=8;
            e->high=(e->high<<8)+255;
        }
    }

    uint32_t bcmeDecodeDirectBits(bcmEncoder *e, int N)
    {
        uint32_t x=0;

        for (int i=0; i<N; ++i)
        {
            const uint32_t mid=e->low+((e->high-e->low)>>1);
            if (e->code<=mid)
            {
                e->high=mid;
                x+=x+1;
            }
            else
            {
                e->low=mid+1;
                x+=x;
            }

            if ((e->low^e->high)<(1<<24))
            {
                e->low<<=8;
                e->high=(e->high<<8)+255;
                e->code=(e->code<<8)+mgetc(g_in);
            }
        }

        return x;
    }

    int bcmeDecodeBit(bcmEncoder *e, uint32_t p)
    {
#if BCM_64BITS
        const uint32_t mid=e->low+(((uint64_t)(e->high-e->low)*p)>>18);
#else
        const uint32_t mid=e->low+(((uint64_t)(e->high-e->low)*(p<<(32-18)))>>32);
#endif
        const int bit=(e->code<=mid);
        if (bit)
            e->high=mid;
        else
            e->low=mid+1;

        while ((e->low^e->high)<(1<<24))
        {
            e->low<<=8;
            e->high=(e->high<<8)+255;
            e->code=(e->code<<8)+mgetc(g_in);
        }

        return bit;
    }

#define BCM_COUNTER_TEMPLATE(RATE) \
typedef struct bcmCounter##RATE { uint16_t p; } bcmCounter##RATE; \
void bcmCounter##RATE##Ctor(bcmCounter##RATE *c) { c->p=1<<15; /* 0.5 */ } \
void bcmCounter##RATE##UpdateBit0(bcmCounter##RATE *c) { c->p-=c->p>>RATE; } \
void bcmCounter##RATE##UpdateBit1(bcmCounter##RATE *c) { c->p+=(c->p^0xFFFF)>>RATE; }

BCM_COUNTER_TEMPLATE(2);
BCM_COUNTER_TEMPLATE(4);
BCM_COUNTER_TEMPLATE(6);

typedef struct bcmCM {
    bcmEncoder enc;
    bcmCounter2 counter0[256];
    bcmCounter4 counter1[256][256];
    bcmCounter6 counter2[2][256][17];
    int c1;
    int c2;
    int run;
} bcmCM;

    void bcmCMCtor(bcmCM *c)
    {
        bcmeCtor(&c->enc);
        for(int i = 0; i < 256; ++i) {
            bcmCounter2Ctor(&c->counter0[i]);
            for(int j = 0; j < 256; ++j) {
                bcmCounter4Ctor(&c->counter1[i][j]);
            }
            for(int k = 0; k < 17; ++k) {
                    bcmCounter6Ctor(&c->counter2[0][i][k]);
                    bcmCounter6Ctor(&c->counter2[1][i][k]);
            }
        }

        c->c1=0;
        c->c2=0;
        c->run=0;

        for (int i=0; i<2; ++i)
        {
            for (int j=0; j<256; ++j)
            {
                for (int k=0; k<17; ++k)
                    c->counter2[i][j][k].p=(k<<12)-(k==16);
            }
        }
    }

    void bcmCMEncode(bcmCM *c, int ch)
    {
        if (c->c1==c->c2)
            ++c->run;
        else
            c->run=0;
        const int f=(c->run>2);

        int ctx=1;
        while (ctx<256)
        {
            const int p0=c->counter0[ctx].p;
            const int p1=c->counter1[c->c1][ctx].p;
            const int p2=c->counter1[c->c2][ctx].p;
            const int p=(((p0+p1)*7)+p2+p2)>>4;

            const int j=p>>12;
            const int x1=c->counter2[f][ctx][j].p;
            const int x2=c->counter2[f][ctx][j+1].p;
            const int ssep=x1+(((x2-x1)*(p&4095))>>12);

            if (ch&128)
            {
                bcmeEncodeBit1(&c->enc, (ssep*3)+p);
                bcmCounter2UpdateBit1(&c->counter0[ctx]);
                bcmCounter4UpdateBit1(&c->counter1[c->c1][ctx]);
                bcmCounter6UpdateBit1(&c->counter2[f][ctx][j]);
                bcmCounter6UpdateBit1(&c->counter2[f][ctx][j+1]);
                ctx+=ctx+1;
            }
            else
            {
                bcmeEncodeBit0(&c->enc, (ssep*3)+p);
                bcmCounter2UpdateBit0(&c->counter0[ctx]);
                bcmCounter4UpdateBit0(&c->counter1[c->c1][ctx]);
                bcmCounter6UpdateBit0(&c->counter2[f][ctx][j]);
                bcmCounter6UpdateBit0(&c->counter2[f][ctx][j+1]);
                ctx+=ctx;
            }

            ch+=ch;
        }

        c->c2=c->c1;
        c->c1=ctx-256;
    }

    int bcmCMDecode(bcmCM *c)
    {
        if (c->c1==c->c2)
            ++c->run;
        else
            c->run=0;
        const int f=(c->run>2);

        int ctx=1;
        while (ctx<256)
        {
            const int p0=c->counter0[ctx].p;
            const int p1=c->counter1[c->c1][ctx].p;
            const int p2=c->counter1[c->c2][ctx].p;
            const int p=(((p0+p1)*7)+p2+p2)>>4;

            const int j=p>>12;
            const int x1=c->counter2[f][ctx][j].p;
            const int x2=c->counter2[f][ctx][j+1].p;
            const int ssep=x1+(((x2-x1)*(p&4095))>>12);

            if (bcmeDecodeBit(&c->enc, (ssep*3)+p))
            {
                bcmCounter2UpdateBit1(&c->counter0[ctx]);
                bcmCounter4UpdateBit1(&c->counter1[c->c1][ctx]);
                bcmCounter6UpdateBit1(&c->counter2[f][ctx][j]);
                bcmCounter6UpdateBit1(&c->counter2[f][ctx][j+1]);
                ctx+=ctx+1;
            }
            else
            {
                bcmCounter2UpdateBit0(&c->counter0[ctx]);
                bcmCounter4UpdateBit0(&c->counter1[c->c1][ctx]);
                bcmCounter6UpdateBit0(&c->counter2[f][ctx][j]);
                bcmCounter6UpdateBit0(&c->counter2[f][ctx][j+1]);
                ctx+=ctx;
            }
        }

        c->c2=c->c1;
        return c->c1=ctx-256;
    }

unsigned bcm_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned level)
{
    mfile infile; minit(&infile, in, inlen); g_in = &infile;
    mfile outfile; minit(&outfile, out, outlen); g_out = &outfile;

    bcmCM cm;
    bcmCMCtor(&cm);

    const int config_tab[10]=
    {
        1<<19,      // -0 - 512KiB, @rlyeh: originally was: 0
        1<<20,      // -1 - 1 MB
        1<<22,      // -2 - 4 MB
        1<<23,      // -3 - 8 MB
        0x00FFFFFF, // -4 - ~16 MB (Default)
        1<<25,      // -5 - 32 MB
        1<<26,      // -6 - 64 MB
        1<<27,      // -7 - 128 MB
        1<<28,      // -8 - 256 MB
        0x7FFFFFFF, // -9 - ~2 GB
    };

    int block_size=config_tab[level];

    int64_t file_size = (int64_t)inlen;
    if (file_size>0 && block_size>file_size)
        block_size=(int)(file_size);

    uint8_t* buf=(uint8_t*)BCM_REALLOC(0, sizeof(uint8_t) * block_size);
    int* ptr=(int*)BCM_REALLOC(0, sizeof(int) * block_size);

    int n;
    while ((n=mread(g_in, buf, block_size))>0)
    {
        const int idx=bcm_divbwt(buf, buf, ptr, n);
        if (idx<1) return 0; // divbwt() failed

        bcmeEncodeDirectBits(&cm.enc, 32, n);
        bcmeEncodeDirectBits(&cm.enc, 32, idx);

        for (int i=0; i<n; ++i)
            bcmCMEncode(&cm, buf[i]);
    }

    bcmeEncodeDirectBits(&cm.enc, 32, 0); // EOF
    // bcmeEncodeDirectBits(&cm.enc, 32, crc32);

    bcmeFlush(&cm.enc);

    BCM_REALLOC(buf, 0); // free
    BCM_REALLOC(ptr, 0); // free

    return mtell(g_out);
}

unsigned bcm_decode(const void *in, unsigned inlen, void *out, unsigned outlen)
{
    mfile infile; minit(&infile, in, inlen); g_in = &infile;
    mfile outfile; minit(&outfile, out, outlen); g_out = &outfile;

    bcmCM cm;
    bcmCMCtor(&cm);

    bcmeInit(&cm.enc);

    int block_size=0;
    uint8_t* buf=NULL;
    uint32_t* ptr=NULL;

    int n;
    while ((n=bcmeDecodeDirectBits(&cm.enc, 32))>0)
    {
        if (block_size==0)
        {
            if ((block_size=n)>=(1<<24)) // 5*N
                buf=(uint8_t*)BCM_REALLOC(0, sizeof(uint8_t) * block_size);

            ptr=(uint32_t*)BCM_REALLOC(0, sizeof(uint32_t) * block_size);
        }

        const int idx=bcmeDecodeDirectBits(&cm.enc, 32);
        if (n>block_size || idx<1 || idx>n) return 0; // corrupt input

        // Inverse BW-transform
        if (n>=(1<<24)) // 5*N
        {
            int t[257]={0};
            for (int i=0; i<n; ++i)
                ++t[(buf[i]=bcmCMDecode(&cm))+1];
            for (int i=1; i<256; ++i)
                t[i]+=t[i-1];
            for (int i=0; i<n; ++i)
                ptr[t[buf[i]]++]=i+(i>=idx);
            for (int p=idx; p;)
            {
                p=ptr[p-1];
                const int c=buf[p-(p>=idx)];
                mputc(g_out, c);
            }
        }
        else // 4*N
        {
            int t[257]={0};
            for (int i=0; i<n; ++i)
                ++t[(ptr[i]=bcmCMDecode(&cm))+1];
            for (int i=1; i<256; ++i)
                t[i]+=t[i-1];
            for (int i=0; i<n; ++i)
                ptr[t[ptr[i]&255]++]|=(i+(i>=idx))<<8;
            for (int p=idx; p;)
            {
                p=ptr[p-1]>>8;
                const int c=ptr[p-(p>=idx)]&255;
                mputc(g_out, c);
            }
        }
    }

    // if (bcmeDecodeDirectBits(&cm.enc, 32)!=crc32) return 0; // crc error

    BCM_REALLOC(buf, 0); // free
    BCM_REALLOC(ptr, 0); // free

    return mtell(g_out);
}

unsigned bcm_bounds(unsigned inlen, unsigned flags) {
        return inlen * 2; // @todo: check src
}

#endif // BCM_C

#line 1 "amalgamated_crush.c" 
// crush.cpp
// Written and placed in the public domain by Ilya Muravyov
// Additional code by @r-lyeh (public domain). @todo: honor unused args inlen/outlen

unsigned crush_encode(const void* in, unsigned inlen, void* out, unsigned outlen, unsigned flags); // [0..(4)..10]
unsigned crush_decode(const void* in, unsigned inlen, void* out, unsigned outlen);
unsigned crush_bounds(unsigned inlen, unsigned flags);


#ifdef CRUSH_C
#pragma once
#ifdef _MSC_VER
#define _CRT_SECECURE_NO_WARNINGS
#define _CRT_DISABLE_PERFCRIT_LOCKS
#endif
#include <stdint.h>
#include <stdlib.h>

// Bit I/O
//
typedef struct bits {
    const uint8_t* g_inbuf;
    uint8_t* g_outbuf;
    int g_inbuf_pos;
    int g_outbuf_pos;
    int bit_buf;
    int bit_count;
} bits;

void bits_init(bits *b, const uint8_t* inbuf, uint8_t* outbuf) {
    b->bit_count=b->bit_buf=b->g_inbuf_pos=b->g_outbuf_pos=0;
    b->g_inbuf = inbuf;
    b->g_outbuf = outbuf;
}

void bits_put(bits *b, int n, int x) {
    b->bit_buf|=x<<b->bit_count;
    b->bit_count+=n;
    while (b->bit_count>=8)
    {
        b->g_outbuf[b->g_outbuf_pos++] = b->bit_buf;
        b->bit_buf>>=8;
        b->bit_count-=8;
    }
}

void bits_flush(bits *b) {
    bits_put(b, 7, 0);
    b->bit_count=b->bit_buf=0;
}

int bits_get(bits *b, int n) {
    while (b->bit_count<n)
    {
        b->bit_buf|=b->g_inbuf[b->g_inbuf_pos++]<<b->bit_count;
        b->bit_count+=8;
    }
    const int x=b->bit_buf&((1<<n)-1);
    b->bit_buf>>=n;
    b->bit_count-=n;
    return x;
}

// LZ77
//

enum { W_BITS=21 }; // Window size (17..23)
enum { W_SIZE=1<<W_BITS };
enum { W_MASK=W_SIZE-1 };
enum { SLOT_BITS=4 };
enum { NUM_SLOTS=1<<SLOT_BITS };

enum { A_BITS=2 }; // 1 xx
enum { B_BITS=2 }; // 01 xx
enum { C_BITS=2 }; // 001 xx
enum { D_BITS=3 }; // 0001 xxx
enum { E_BITS=5 }; // 00001 xxxxx
enum { F_BITS=9 }; // 00000 xxxxxxxxx
enum { A=1<<A_BITS };
enum { B=(1<<B_BITS)+A };
enum { C=(1<<C_BITS)+B };
enum { D=(1<<D_BITS)+C };
enum { E=(1<<E_BITS)+D };
enum { F=(1<<F_BITS)+E };
enum { MIN_MATCH=3 };
enum { MAX_MATCH=(F-1)+MIN_MATCH };

enum { TOO_FAR=1<<16 };

enum { HASH1_LEN=MIN_MATCH };
enum { HASH2_LEN=MIN_MATCH+1 };
enum { HASH1_BITS=21 };
enum { HASH2_BITS=24 };
enum { HASH1_SIZE=1<<HASH1_BITS };
enum { HASH2_SIZE=1<<HASH2_BITS };
enum { HASH1_MASK=HASH1_SIZE-1 };
enum { HASH2_MASK=HASH2_SIZE-1 };
enum { HASH1_SHIFT=(HASH1_BITS+(HASH1_LEN-1))/HASH1_LEN };
enum { HASH2_SHIFT=(HASH2_BITS+(HASH2_LEN-1))/HASH2_LEN };

static inline int update_hash1(int h, int c) {
    return ((h<<HASH1_SHIFT)+c)&HASH1_MASK;
}

static inline int update_hash2(int h, int c) {
    return ((h<<HASH2_SHIFT)+c)&HASH2_MASK;
}

static inline int get_min(int a, int b) {
    return a<b?a:b;
}

static inline int get_max(int a, int b) {
    return a>b?a:b;
}

static inline int get_penalty(int a, int b) {
    int p=0;
    while (a>b)
    {
        a>>=3;
        ++p;
    }
    return p;
}

static size_t crush_compress(const uint8_t* buf, size_t size, uint8_t* outbuf, size_t outlen, size_t level) {
    static int head[HASH1_SIZE+HASH2_SIZE];
    static int prev[W_SIZE];

    //const int max_chain[]={4, 256, 1<<12}; // original [0fast..2uber]
    const int max_chain[11] = { 0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 1<<12 }; //[0fastest..10uber]
    const int max_level = sizeof(max_chain)/sizeof(max_chain[0]);
    level = level > max_level ? max_level : level;

    bits bits;

    {
        for (int i=0; i<HASH1_SIZE+HASH2_SIZE; ++i)
            head[i]=-1;

        int h1=0;
        int h2=0;
        for (int i=0; i<HASH1_LEN; ++i)
            h1=update_hash1(h1, buf[i]);
        for (int i=0; i<HASH2_LEN; ++i)
            h2=update_hash2(h2, buf[i]);

        bits_init(&bits, NULL, outbuf);

        size_t p=0;
        while (p<size)
        {
            int len=MIN_MATCH-1;
            int offset=W_SIZE;

            const int max_match=get_min(MAX_MATCH, size-p);
            const int limit=get_max(p-W_SIZE, 0);

            if (head[h1]>=limit)
            {
                int s=head[h1];
                if (buf[s]==buf[p])
                {
                    int l=0;
                    while (++l<max_match)
                        if (buf[s+l]!=buf[p+l])
                            break;
                    if (l>len)
                    {
                        len=l;
                        offset=p-s;
                    }
                }
            }

            if (len<MAX_MATCH)
            {
                int chain_len=max_chain[level];
                int s=head[h2+HASH1_SIZE];

                while ((chain_len--!=0)&&(s>=limit))
                {
                    if ((buf[s+len]==buf[p+len])&&(buf[s]==buf[p]))
                    {   
                        int l=0;
                        while (++l<max_match)
                            if (buf[s+l]!=buf[p+l])
                                break;
                        if (l>len+get_penalty((p-s)>>4, offset))
                        {
                            len=l;
                            offset=p-s;
                        }
                        if (l==max_match)
                            break;
                    }
                    s=prev[s&W_MASK];
                }
            }

            if ((len==MIN_MATCH)&&(offset>TOO_FAR))
                len=0;

            if ((level>=2)&&(len>=MIN_MATCH)&&(len<max_match))
            {
                const int next_p=p+1;
                const int max_lazy=get_min(len+4, max_match);

                int chain_len=max_chain[level];
                int s=head[update_hash2(h2, buf[next_p+(HASH2_LEN-1)])+HASH1_SIZE];

                while ((chain_len--!=0)&&(s>=limit))
                {
                    if ((buf[s+len]==buf[next_p+len])&&(buf[s]==buf[next_p]))
                    {
                        int l=0;
                        while (++l<max_lazy)
                            if (buf[s+l]!=buf[next_p+l])
                                break;
                        if (l>len+get_penalty(next_p-s, offset))
                        {
                            len=0;
                            break;
                        }
                        if (l==max_lazy)
                            break;
                    }
                    s=prev[s&W_MASK];
                }
            }

            if (len>=MIN_MATCH) // Match
            {
                bits_put(&bits, 1, 1);

                const int l=len-MIN_MATCH;
                if (l<A)
                {
                    bits_put(&bits, 1, 1); // 1
                    bits_put(&bits, A_BITS, l);
                }
                else if (l<B)
                {
                    bits_put(&bits, 2, 1<<1); // 01
                    bits_put(&bits, B_BITS, l-A);
                }
                else if (l<C)
                {
                    bits_put(&bits, 3, 1<<2); // 001
                    bits_put(&bits, C_BITS, l-B);
                }
                else if (l<D)
                {
                    bits_put(&bits, 4, 1<<3); // 0001
                    bits_put(&bits, D_BITS, l-C);
                }
                else if (l<E)
                {
                    bits_put(&bits, 5, 1<<4); // 00001
                    bits_put(&bits, E_BITS, l-D);
                }
                else
                {
                    bits_put(&bits, 5, 0); // 00000
                    bits_put(&bits, F_BITS, l-E);
                }

                --offset;
                int log=W_BITS-NUM_SLOTS;
                while (offset>=(2<<log))
                    ++log;
                bits_put(&bits, SLOT_BITS, log-(W_BITS-NUM_SLOTS));
                if (log>(W_BITS-NUM_SLOTS))
                    bits_put(&bits, log, offset-(1<<log));
                else
                    bits_put(&bits, W_BITS-(NUM_SLOTS-1), offset);
            }
            else // Literal
            {
                len=1;
                bits_put(&bits, 9, buf[p]<<1); // 0 xxxxxxxx
            }

            while (len--!=0) // Insert new strings
            {
                head[h1]=p;
                prev[p&W_MASK]=head[h2+HASH1_SIZE];
                head[h2+HASH1_SIZE]=p;
                ++p;
                h1=update_hash1(h1, buf[p+(HASH1_LEN-1)]);
                h2=update_hash2(h2, buf[p+(HASH2_LEN-1)]);
            }
        }

        bits_flush(&bits);
    }

    return bits.g_outbuf_pos;
}

static size_t crush_decompress(const uint8_t* inbuf, size_t inlen, uint8_t* outbuf, size_t outsize)
{
    if (inlen<1)
    {
        //fprintf(stderr, "Corrupted stream: size=%d\n", (int)inlen);
        return 0;
    }

    bits bits;
    bits_init(&bits, inbuf, NULL);

    int p=0;
    while (bits.g_inbuf_pos<(int)inlen)
    {
        if (bits_get(&bits, 1))
        {
            int len;
            /**/ if (bits_get(&bits, 1)) len=bits_get(&bits, A_BITS);
            else if (bits_get(&bits, 1)) len=bits_get(&bits, B_BITS)+A;
            else if (bits_get(&bits, 1)) len=bits_get(&bits, C_BITS)+B;
            else if (bits_get(&bits, 1)) len=bits_get(&bits, D_BITS)+C;
            else if (bits_get(&bits, 1)) len=bits_get(&bits, E_BITS)+D;
            else                         len=bits_get(&bits, F_BITS)+E;

            const int log=bits_get(&bits, SLOT_BITS)+(W_BITS-NUM_SLOTS);
            int s=~(log>(W_BITS-NUM_SLOTS)
                ?bits_get(&bits, log)+(1<<log)
                :bits_get(&bits, W_BITS-(NUM_SLOTS-1)))+p;
            if (s<0)
            {
                //fprintf(stderr, "Corrupted stream: s=%d p=%d inlen=%d\n", s, p, (int)inlen);
                return 0;
            }

            outbuf[p++]=outbuf[s++];
            outbuf[p++]=outbuf[s++];
            outbuf[p++]=outbuf[s++];
            while (len--!=0)
                outbuf[p++]=outbuf[s++];
        }
        else
            outbuf[p++]=bits_get(&bits, 8);
    }

    return p;
}


unsigned crush_encode(const void* in, unsigned inlen, void* out, unsigned outlen, unsigned flags) {
    unsigned level = flags > 10 ? 10 : flags < 0 ? 0 : flags;
    return crush_compress((const uint8_t*)in, (size_t)inlen, (uint8_t*)out, (size_t)outlen, (size_t)level);
}
unsigned crush_decode(const void* in, unsigned inlen, void* out, unsigned outlen) {
    return crush_decompress((const uint8_t*)in, (size_t)inlen, (uint8_t*)out, (size_t)outlen);
}
unsigned crush_bounds(unsigned inlen, unsigned flags) {
    return (unsigned)(inlen * 1.1) + 16; // @todo: check src
}

#endif // CRUSH_C

#ifdef CRUSH_DEMO
#pragma once
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level = 1;
    char out[128];
    size_t outlen = crush_encode(longcopy, strlen(longcopy)+1, out, 128, level );
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    size_t unpacked = crush_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // CRUSH_DEMO

#line 1 "amalgamated_deflate.c" 
// miniz.c v1.15 r4 - public domain de/inflate. See "unlicense" statement at http://unlicense.org/
// Rich Geldreich <richgel99@gmail.com>, last updated Oct. 13, 2013. Then stripped down by @r-lyeh.
// Implements RFC 1950: http://www.ietf.org/rfc/rfc1950.txt and RFC 1951: http://www.ietf.org/rfc/rfc1951.txt

unsigned deflate_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags); // [0..(6)..9][10 (uber)]
unsigned deflate_decode(const void *in, unsigned inlen, void *out, unsigned outlen);
unsigned deflate_bounds(unsigned inlen, unsigned flags);


#ifdef DEFLATE_C
#pragma once

#define SDEFL_IMPLEMENTATION
#define SINFL_IMPLEMENTATION

/*
# Small Deflate
`sdefl` is a small bare bone lossless compression library in ANSI C (ISO C90)
which implements the Deflate (RFC 1951) compressed data format specification standard.
It is mainly tuned to get as much speed and compression ratio from as little code
as needed to keep the implementation as concise as possible.

## Features
- Portable single header and source file duo written in ANSI C (ISO C90)
- Dual license with either MIT or public domain
- Small implementation
    - Deflate: 525 LoC
    - Inflate: 320 LoC
- Webassembly:
    - Deflate ~3.7 KB (~2.2KB compressed)
    - Inflate ~3.6 KB (~2.2KB compressed)

## Usage:
This file behaves differently depending on what symbols you define
before including it.

Header-File mode:
If you do not define `SDEFL_IMPLEMENTATION` before including this file, it
will operate in header only mode. In this mode it declares all used structs
and the API of the library without including the implementation of the library.

Implementation mode:
If you define `SDEFL_IMPLEMENTATION` before including this file, it will
compile the implementation . Make sure that you only include
this file implementation in *one* C or C++ file to prevent collisions.

### Benchmark

| Compressor name         | Compression| Decompress.| Compr. size | Ratio |
| ------------------------| -----------| -----------| ----------- | ----- |
| sdefl 1.0 -0            |   127 MB/s |   233 MB/s |    40004116 | 39.88 |
| sdefl 1.0 -1            |   111 MB/s |   259 MB/s |    38940674 | 38.82 |
| sdefl 1.0 -5            |    45 MB/s |   275 MB/s |    36577183 | 36.46 |
| sdefl 1.0 -7            |    38 MB/s |   276 MB/s |    36523781 | 36.41 |
| zlib 1.2.11 -1          |    72 MB/s |   307 MB/s |    42298774 | 42.30 |
| zlib 1.2.11 -6          |    24 MB/s |   313 MB/s |    36548921 | 36.55 |
| zlib 1.2.11 -9          |    20 MB/s |   314 MB/s |    36475792 | 36.48 |
| miniz 1.0 -1            |   122 MB/s |   208 MB/s |    48510028 | 48.51 |
| miniz 1.0 -6            |    27 MB/s |   260 MB/s |    36513697 | 36.51 |
| miniz 1.0 -9            |    23 MB/s |   261 MB/s |    36460101 | 36.46 |
| libdeflate 1.3 -1       |   147 MB/s |   667 MB/s |    39597378 | 39.60 |
| libdeflate 1.3 -6       |    69 MB/s |   689 MB/s |    36648318 | 36.65 |
| libdeflate 1.3 -9       |    13 MB/s |   672 MB/s |    35197141 | 35.20 |
| libdeflate 1.3 -12      |  8.13 MB/s |   670 MB/s |    35100568 | 35.10 |

### Compression
Results on the [Silesia compression corpus](http://sun.aei.polsl.pl/~sdeor/index.php?page=silesia):

| File    |   Original | `sdefl 0`      | `sdefl 5`     | `sdefl 7` |
| :------ | ---------: | -----------------: | ---------: | ----------: |
| dickens | 10.192.446 |  4,260,187|  3,845,261|   3,833,657 |
| mozilla | 51.220.480 | 20,774,706 | 19,607,009 |  19,565,867 |
| mr      |  9.970.564 | 3,860,531 |  3,673,460 |   3,665,627 |
| nci     | 33.553.445 | 4,030,283 |  3,094,526 |   3,006,075 |
| ooffice |  6.152.192 | 3,320,063 |  3,186,373 |   3,183,815 |
| osdb    | 10.085.684 | 3,919,646 |  3,649,510 |   3,649,477 |
| reymont |  6.627.202 | 2,263,378 |  1,857,588 |   1,827,237 |
| samba   | 21.606.400 | 6,121,797 |  5,462,670 |   5,450,762 |
| sao     |  7.251.944 | 5,612,421 |  5,485,380 |   5,481,765 |
| webster | 41.458.703 | 13,972,648 | 12,059,432 |  11,991,421 |
| xml     |  5.345.280 | 886,620|    674,009 |     662,141 |
| x-ray   |  8.474.240 | 6,304,655 |  6,244,779 |   6,244,779 |

## License
```
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2020 Micha Mettke
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
```
*/
#ifndef SDEFL_H_INCLUDED
#define SDEFL_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#define SDEFL_MAX_OFF   (1 << 15)
#define SDEFL_WIN_SIZ   SDEFL_MAX_OFF
#define SDEFL_WIN_MSK   (SDEFL_WIN_SIZ-1)

#define SDEFL_HASH_BITS 15
#define SDEFL_HASH_SIZ  (1 << SDEFL_HASH_BITS)
#define SDEFL_HASH_MSK  (SDEFL_HASH_SIZ-1)

#define SDEFL_MIN_MATCH 4
#define SDEFL_BLK_MAX   (256*1024)
#define SDEFL_SEQ_SIZ   ((SDEFL_BLK_MAX + SDEFL_MIN_MATCH)/SDEFL_MIN_MATCH)

#define SDEFL_SYM_MAX   (288)
#define SDEFL_OFF_MAX   (32)
#define SDEFL_PRE_MAX   (19)

#define SDEFL_LVL_MIN   0
#define SDEFL_LVL_DEF   5
#define SDEFL_LVL_MAX   8

struct sdefl_freq {
  unsigned lit[SDEFL_SYM_MAX];
  unsigned off[SDEFL_OFF_MAX];
};
struct sdefl_code_words {
  unsigned lit[SDEFL_SYM_MAX];
  unsigned off[SDEFL_OFF_MAX];
};
struct sdefl_lens {
  unsigned char lit[SDEFL_SYM_MAX];
  unsigned char off[SDEFL_OFF_MAX];
};
struct sdefl_codes {
  struct sdefl_code_words word;
  struct sdefl_lens len;
};
struct sdefl_seqt {
  int off, len;
};
struct sdefl {
  int bits, bitcnt;
  int tbl[SDEFL_HASH_SIZ];
  int prv[SDEFL_WIN_SIZ];

  int seq_cnt;
  struct sdefl_seqt seq[SDEFL_SEQ_SIZ];
  struct sdefl_freq freq;
  struct sdefl_codes cod;
};
extern int sdefl_bound(int in_len);
extern int sdeflate(struct sdefl *s, void *o, const void *i, int n, int lvl);
extern int zsdeflate(struct sdefl *s, void *o, const void *i, int n, int lvl);

#ifdef __cplusplus
}
#endif

#endif /* SDEFL_H_INCLUDED */

#ifdef SDEFL_IMPLEMENTATION

#include <assert.h> /* assert */
#include <string.h> /* memcpy */
#include <limits.h> /* CHAR_BIT */

#define SDEFL_NIL               (-1)
#define SDEFL_MAX_MATCH         258
#define SDEFL_MAX_CODE_LEN      (15)
#define SDEFL_SYM_BITS          (10u)
#define SDEFL_SYM_MSK           ((1u << SDEFL_SYM_BITS)-1u)
#define SDEFL_LIT_LEN_CODES     (14)
#define SDEFL_OFF_CODES         (15)
#define SDEFL_PRE_CODES         (7)
#define SDEFL_CNT_NUM(n)        ((((n)+3u/4u)+3u)&~3u)
#define SDEFL_EOB               (256)

#define sdefl_npow2(n) (1 << (sdefl_ilog2((n)-1) + 1))

static int
sdefl_ilog2(int n) {
  if (!n) return 0;
#ifdef _MSC_VER
  unsigned long msbp = 0;
  _BitScanReverse(&msbp, (unsigned long)n);
  return (int)msbp;
#elif defined(__GNUC__) || defined(__clang__)
  return (int)sizeof(unsigned long) * CHAR_BIT - 1 - __builtin_clzl((unsigned long)n);
#else
  #define lt(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
  static const char tbl[256] = {
    0,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3,lt(4), lt(5), lt(5), lt(6), lt(6), lt(6), lt(6),
    lt(7), lt(7), lt(7), lt(7), lt(7), lt(7), lt(7), lt(7)};
  int tt, t;
  if ((tt = (n >> 16))) {
    return (t = (tt >> 8)) ? 24 + tbl[t] : 16 + tbl[tt];
  } else {
    return (t = (n >> 8)) ? 8 + tbl[t] : tbl[n];
  }
  #undef lt
#endif
}
static unsigned
sdefl_uload32(const void *p) {
  /* hopefully will be optimized to an unaligned read */
  unsigned n = 0;
  memcpy(&n, p, sizeof(n));
  return n;
}
static unsigned
sdefl_hash32(const void *p) {
  unsigned n = sdefl_uload32(p);
  return (n * 0x9E377989) >> (32 - SDEFL_HASH_BITS);
}
static void
sdefl_put(unsigned char **dst, struct sdefl *s, int code, int bitcnt) {
  s->bits |= (code << s->bitcnt);
  s->bitcnt += bitcnt;
  while (s->bitcnt >= 8) {
    unsigned char *tar = *dst;
    *tar = (unsigned char)(s->bits & 0xFF);
    s->bits >>= 8;
    s->bitcnt -= 8;
    *dst = *dst + 1;
  }
}
static void
sdefl_heap_sub(unsigned A[], unsigned len, unsigned sub) {
  unsigned c, p = sub;
  unsigned v = A[sub];
  while ((c = p << 1) <= len) {
    if (c < len && A[c + 1] > A[c]) c++;
    if (v >= A[c]) break;
    A[p] = A[c], p = c;
  }
  A[p] = v;
}
static void
sdefl_heap_array(unsigned *A, unsigned len) {
  unsigned sub;
  for (sub = len >> 1; sub >= 1; sub--)
    sdefl_heap_sub(A, len, sub);
}
static void
sdefl_heap_sort(unsigned *A, unsigned n) {
  A--;
  sdefl_heap_array(A, n);
  while (n >= 2) {
    unsigned tmp = A[n];
    A[n--] = A[1];
    A[1] = tmp;
    sdefl_heap_sub(A, n, 1);
  }
}
static unsigned
sdefl_sort_sym(unsigned sym_cnt, unsigned *freqs,
               unsigned char *lens, unsigned *sym_out) {
  unsigned cnts[SDEFL_CNT_NUM(SDEFL_SYM_MAX)] = {0};
  unsigned cnt_num = SDEFL_CNT_NUM(sym_cnt);
  unsigned used_sym = 0;
  unsigned sym, i;
  for (sym = 0; sym < sym_cnt; sym++)
    cnts[freqs[sym] < cnt_num-1 ? freqs[sym]: cnt_num-1]++;
  for (i = 1; i < cnt_num; i++) {
    unsigned cnt = cnts[i];
    cnts[i] = used_sym;
    used_sym += cnt;
  }
  for (sym = 0; sym < sym_cnt; sym++) {
    unsigned freq = freqs[sym];
    if (freq) {
        unsigned idx = freq < cnt_num-1 ? freq : cnt_num-1;
        sym_out[cnts[idx]++] = sym | (freq << SDEFL_SYM_BITS);
    } else lens[sym] = 0;
  }
  sdefl_heap_sort(sym_out + cnts[cnt_num-2], cnts[cnt_num-1] - cnts[cnt_num-2]);
  return used_sym;
}
static void
sdefl_build_tree(unsigned *A, unsigned sym_cnt) {
  unsigned i = 0, b = 0, e = 0;
  do {
    unsigned m, n, freq_shift;
    if (i != sym_cnt && (b == e || (A[i] >> SDEFL_SYM_BITS) <= (A[b] >> SDEFL_SYM_BITS)))
      m = i++;
    else m = b++;
    if (i != sym_cnt && (b == e || (A[i] >> SDEFL_SYM_BITS) <= (A[b] >> SDEFL_SYM_BITS)))
      n = i++;
    else n = b++;

    freq_shift = (A[m] & ~SDEFL_SYM_MSK) + (A[n] & ~SDEFL_SYM_MSK);
    A[m] = (A[m] & SDEFL_SYM_MSK) | (e << SDEFL_SYM_BITS);
    A[n] = (A[n] & SDEFL_SYM_MSK) | (e << SDEFL_SYM_BITS);
    A[e] = (A[e] & SDEFL_SYM_MSK) | freq_shift;
  } while (sym_cnt - ++e > 1);
}
static void
sdefl_gen_len_cnt(unsigned *A, unsigned root, unsigned *len_cnt,
                  unsigned max_code_len) {
  int n;
  unsigned i;
  for (i = 0; i <= max_code_len; i++)
    len_cnt[i] = 0;
  len_cnt[1] = 2;

  A[root] &= SDEFL_SYM_MSK;
  for (n = (int)root - 1; n >= 0; n--) {
    unsigned p = A[n] >> SDEFL_SYM_BITS;
    unsigned pdepth = A[p] >> SDEFL_SYM_BITS;
    unsigned depth = pdepth + 1;
    unsigned len = depth;

    A[n] = (A[n] & SDEFL_SYM_MSK) | (depth << SDEFL_SYM_BITS);
    if (len >= max_code_len) {
      len = max_code_len;
      do len--; while (!len_cnt[len]);
    }
    len_cnt[len]--;
    len_cnt[len+1] += 2;
  }
}
static void
sdefl_gen_codes(unsigned *A, unsigned char *lens, const unsigned *len_cnt,
                unsigned max_code_word_len, unsigned sym_cnt) {
  unsigned i, sym, len, nxt[SDEFL_MAX_CODE_LEN + 1];
  for (i = 0, len = max_code_word_len; len >= 1; len--) {
    unsigned cnt = len_cnt[len];
    while (cnt--) lens[A[i++] & SDEFL_SYM_MSK] = (unsigned char)len;
  }
  nxt[0] = nxt[1] = 0;
  for (len = 2; len <= max_code_word_len; len++)
    nxt[len] = (nxt[len-1] + len_cnt[len-1]) << 1;
  for (sym = 0; sym < sym_cnt; sym++)
    A[sym] = nxt[lens[sym]]++;
}
static unsigned
sdefl_rev(unsigned c, unsigned char n) {
  c = ((c & 0x5555) << 1) | ((c & 0xAAAA) >> 1);
  c = ((c & 0x3333) << 2) | ((c & 0xCCCC) >> 2);
  c = ((c & 0x0F0F) << 4) | ((c & 0xF0F0) >> 4);
  c = ((c & 0x00FF) << 8) | ((c & 0xFF00) >> 8);
  return c >> (16-n);
}
static void
sdefl_huff(unsigned char *lens, unsigned *codes, unsigned *freqs,
           unsigned num_syms, unsigned max_code_len) {
  unsigned c, *A = codes;
  unsigned len_cnt[SDEFL_MAX_CODE_LEN + 1];
  unsigned used_syms = sdefl_sort_sym(num_syms, freqs, lens, A);
  if (!used_syms) return;
  if (used_syms == 1) {
    unsigned s = A[0] & SDEFL_SYM_MSK;
    unsigned i = s ? s : 1;
    codes[0] = 0, lens[0] = 1;
    codes[i] = 1, lens[i] = 1;
    return;
  }
  sdefl_build_tree(A, used_syms);
  sdefl_gen_len_cnt(A, used_syms-2, len_cnt, max_code_len);
  sdefl_gen_codes(A, lens, len_cnt, max_code_len, num_syms);
  for (c = 0; c < num_syms; c++) {
    codes[c] = sdefl_rev(codes[c], lens[c]);
  }
}
struct sdefl_symcnt {
  int items;
  int lit;
  int off;
};
static void
sdefl_precode(struct sdefl_symcnt *cnt, unsigned *freqs, unsigned *items,
              const unsigned char *litlen, const unsigned char *offlen) {
  unsigned *at = items;
  unsigned run_start = 0;

  unsigned total = 0;
  unsigned char lens[SDEFL_SYM_MAX + SDEFL_OFF_MAX];
  for (cnt->lit = SDEFL_SYM_MAX; cnt->lit > 257; cnt->lit--)
    if (litlen[cnt->lit - 1]) break;
  for (cnt->off = SDEFL_OFF_MAX; cnt->off > 1; cnt->off--)
    if (offlen[cnt->off - 1]) break;

  total = (unsigned)(cnt->lit + cnt->off);
  memcpy(lens, litlen, sizeof(unsigned char) * cnt->lit);
  memcpy(lens + cnt->lit, offlen, sizeof(unsigned char) * cnt->off);
  do {
    unsigned len = lens[run_start];
    unsigned run_end = run_start;
    do run_end++; while (run_end != total && len == lens[run_end]);
    if (!len) {
      while ((run_end - run_start) >= 11) {
        unsigned n = (run_end - run_start) - 11;
        unsigned xbits =  n < 0x7f ? n : 0x7f;
        freqs[18]++;
        *at++ = 18u | (xbits << 5u);
        run_start += 11 + xbits;
      }
      if ((run_end - run_start) >= 3) {
        unsigned n = (run_end - run_start) - 3;
        unsigned xbits =  n < 0x7 ? n : 0x7;
        freqs[17]++;
        *at++ = 17u | (xbits << 5u);
        run_start += 3 + xbits;
      }
    } else if ((run_end - run_start) >= 4) {
      freqs[len]++;
      *at++ = len;
      run_start++;
      do {
        unsigned xbits = (run_end - run_start) - 3;
        xbits = xbits < 0x03 ? xbits : 0x03;
        *at++ = 16 | (xbits << 5);
        run_start += 3 + xbits;
        freqs[16]++;
      } while ((run_end - run_start) >= 3);
    }
    while (run_start != run_end) {
      freqs[len]++;
      *at++ = len;
      run_start++;
    }
  } while (run_start != total);
  cnt->items = (int)(at - items);
}
struct sdefl_match_codes {
  int ls, lc;
  int dc, dx;
};
static void
sdefl_match_codes(struct sdefl_match_codes *cod, int dist, int len) {
  static const short dxmax[] = {0,6,12,24,48,96,192,384,768,1536,3072,6144,12288,24576};
  static const unsigned char lslot[258+1] = {
    0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 12,
    12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16, 16,
    16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18,
    18, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
    21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22,
    22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25,
    25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
    25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
    27, 27, 28
  };
  cod->ls = lslot[len];
  cod->lc = 257 + cod->ls;
  cod->dx = sdefl_ilog2(sdefl_npow2(dist) >> 2);
  cod->dc = cod->dx ? ((cod->dx + 1) << 1) + (dist > dxmax[cod->dx]) : dist-1;
}
static void
sdefl_match(unsigned char **dst, struct sdefl *s, int dist, int len) {
  static const char lxn[] = {0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0};
  static const short lmin[] = {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,
      51,59,67,83,99,115,131,163,195,227,258};
  static const short dmin[] = {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,
      385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};

  struct sdefl_match_codes cod;
  sdefl_match_codes(&cod, dist, len);
  sdefl_put(dst, s, (int)s->cod.word.lit[cod.lc], s->cod.len.lit[cod.lc]);
  sdefl_put(dst, s, len - lmin[cod.ls], lxn[cod.ls]);
  sdefl_put(dst, s, (int)s->cod.word.off[cod.dc], s->cod.len.off[cod.dc]);
  sdefl_put(dst, s, dist - dmin[cod.dc], cod.dx);
}
static void
sdefl_flush(unsigned char **dst, struct sdefl *s, int is_last,
            const unsigned char *in) {
  int j, i = 0, item_cnt = 0;
  struct sdefl_symcnt symcnt = {0};
  unsigned codes[SDEFL_PRE_MAX];
  unsigned char lens[SDEFL_PRE_MAX];
  unsigned freqs[SDEFL_PRE_MAX] = {0};
  unsigned items[SDEFL_SYM_MAX + SDEFL_OFF_MAX];
  static const unsigned char perm[SDEFL_PRE_MAX] = {16,17,18,0,8,7,9,6,10,5,11,
      4,12,3,13,2,14,1,15};

  /* huffman codes */
  s->freq.lit[SDEFL_EOB]++;
  sdefl_huff(s->cod.len.lit, s->cod.word.lit, s->freq.lit, SDEFL_SYM_MAX, SDEFL_LIT_LEN_CODES);
  sdefl_huff(s->cod.len.off, s->cod.word.off, s->freq.off, SDEFL_OFF_MAX, SDEFL_OFF_CODES);
  sdefl_precode(&symcnt, freqs, items, s->cod.len.lit, s->cod.len.off);
  sdefl_huff(lens, codes, freqs, SDEFL_PRE_MAX, SDEFL_PRE_CODES);
  for (item_cnt = SDEFL_PRE_MAX; item_cnt > 4; item_cnt--) {
    if (lens[perm[item_cnt - 1]]) break;
  }
  /* block header */
  sdefl_put(dst, s, is_last ? 0x01 : 0x00, 1); /* block */
  sdefl_put(dst, s, 0x02, 2); /* dynamic huffman */
  sdefl_put(dst, s, symcnt.lit - 257, 5);
  sdefl_put(dst, s, symcnt.off - 1, 5);
  sdefl_put(dst, s, item_cnt - 4, 4);
  for (i = 0; i < item_cnt; ++i)
    sdefl_put(dst, s, lens[perm[i]], 3);
  for (i = 0; i < symcnt.items; ++i) {
    unsigned sym = items[i] & 0x1F;
    sdefl_put(dst, s, (int)codes[sym], lens[sym]);
    if (sym < 16) continue;
    if (sym == 16) sdefl_put(dst, s, items[i] >> 5, 2);
    else if(sym == 17) sdefl_put(dst, s, items[i] >> 5, 3);
    else sdefl_put(dst, s, items[i] >> 5, 7);
  }
  /* block sequences */
  for (i = 0; i < s->seq_cnt; ++i) {
    if (s->seq[i].off >= 0)
      for (j = 0; j < s->seq[i].len; ++j) {
        int c = in[s->seq[i].off + j];
        sdefl_put(dst, s, (int)s->cod.word.lit[c], s->cod.len.lit[c]);
      }
    else sdefl_match(dst, s, -s->seq[i].off, s->seq[i].len);
  }
  sdefl_put(dst, s, (int)(s)->cod.word.lit[SDEFL_EOB], (s)->cod.len.lit[SDEFL_EOB]);
  memset(&s->freq, 0, sizeof(s->freq));
  s->seq_cnt = 0;
}
static void
sdefl_seq(struct sdefl *s, int off, int len) {
  assert(s->seq_cnt + 2 < SDEFL_SEQ_SIZ);
  s->seq[s->seq_cnt].off = off;
  s->seq[s->seq_cnt].len = len;
  s->seq_cnt++;
}
static void
sdefl_reg_match(struct sdefl *s, int off, int len) {
  struct sdefl_match_codes cod;
  sdefl_match_codes(&cod, off, len);
  s->freq.lit[cod.lc]++;
  s->freq.off[cod.dc]++;
}
struct sdefl_match {
  int off;
  int len;
};
static void
sdefl_fnd(struct sdefl_match *m, const struct sdefl *s,
          int chain_len, int max_match, const unsigned char *in, int p) {
  int i = s->tbl[sdefl_hash32(&in[p])];
  int limit = ((p-SDEFL_WIN_SIZ)<SDEFL_NIL)?SDEFL_NIL:(p-SDEFL_WIN_SIZ);
  while (i > limit) {
    if (in[i+m->len] == in[p+m->len] &&
        (sdefl_uload32(&in[i]) == sdefl_uload32(&in[p]))){
      int n = SDEFL_MIN_MATCH;
      while (n < max_match && in[i+n] == in[p+n]) n++;
      if (n > m->len) {
        m->len = n, m->off = p - i;
        if (n == max_match) break;
      }
    }
    if (!(--chain_len)) break;
    i = s->prv[i&SDEFL_WIN_MSK];
  }
}
static int
sdefl_compr(struct sdefl *s, unsigned char *out, const unsigned char *in,
            int in_len, int lvl) {
  unsigned char *q = out;
  static const unsigned char pref[] = {8,10,14,24,30,48,65,96,130};
  int max_chain = (lvl < 8) ? (1 << (lvl + 1)): (1 << 13);
  int n, i = 0, litlen = 0;
  for (n = 0; n < SDEFL_HASH_SIZ; ++n) {
    s->tbl[n] = SDEFL_NIL;
  }
  do {int blk_end = i + SDEFL_BLK_MAX < in_len ? i + SDEFL_BLK_MAX : in_len;
    while (i < blk_end) {
      struct sdefl_match m = {0};
      int max_match = ((in_len-i)>SDEFL_MAX_MATCH) ? SDEFL_MAX_MATCH:(in_len-i);
      int nice_match = pref[lvl] < max_match ? pref[lvl] : max_match;
      int run = 1, inc = 1, run_inc;
      if (max_match > SDEFL_MIN_MATCH) {
        sdefl_fnd(&m, s, max_chain, max_match, in, i);
      }
      if (lvl >= 5 && m.len >= SDEFL_MIN_MATCH && m.len < nice_match){
        struct sdefl_match m2 = {0};
        sdefl_fnd(&m2, s, max_chain, m.len+1, in, i+1);
        m.len = (m2.len > m.len) ? 0 : m.len;
      }
      if (m.len >= SDEFL_MIN_MATCH) {
        if (litlen) {
          sdefl_seq(s, i - litlen, litlen);
          litlen = 0;
        }
        sdefl_seq(s, -m.off, m.len);
        sdefl_reg_match(s, m.off, m.len);
        if (lvl < 2 && m.len >= nice_match) {
          inc = m.len;
        } else {
          run = m.len;
        }
      } else {
        s->freq.lit[in[i]]++;
        litlen++;
      }
      run_inc = run * inc;
      if (in_len - (i + run_inc) > SDEFL_MIN_MATCH) {
        while (run-- > 0) {
          unsigned h = sdefl_hash32(&in[i]);
          s->prv[i&SDEFL_WIN_MSK] = s->tbl[h];
          s->tbl[h] = i, i += inc;
        }
      } else {
        i += run_inc;
      }
    }
    if (litlen) {
      sdefl_seq(s, i - litlen, litlen);
      litlen = 0;
    }
    sdefl_flush(&q, s, blk_end == in_len, in);
  } while (i < in_len);

  if (s->bitcnt)
    sdefl_put(&q, s, 0x00, 8 - s->bitcnt);
  return (int)(q - out);
}
extern int
sdeflate(struct sdefl *s, void *out, const void *in, int n, int lvl) {
  s->bits = s->bitcnt = 0;
  return sdefl_compr(s, (unsigned char*)out, (const unsigned char*)in, n, lvl);
}
static unsigned
sdefl_adler32(unsigned adler32, const unsigned char *in, int in_len) {
  #define SDEFL_ADLER_INIT (1)
  const unsigned ADLER_MOD = 65521;
  unsigned s1 = adler32 & 0xffff;
  unsigned s2 = adler32 >> 16;
  unsigned blk_len, i;

  blk_len = in_len % 5552;
  while (in_len) {
    for (i = 0; i + 7 < blk_len; i += 8) {
      s1 += in[0]; s2 += s1;
      s1 += in[1]; s2 += s1;
      s1 += in[2]; s2 += s1;
      s1 += in[3]; s2 += s1;
      s1 += in[4]; s2 += s1;
      s1 += in[5]; s2 += s1;
      s1 += in[6]; s2 += s1;
      s1 += in[7]; s2 += s1;
      in += 8;
    }
    for (; i < blk_len; ++i) {
      s1 += *in++, s2 += s1;
    }
    s1 %= ADLER_MOD;
    s2 %= ADLER_MOD;
    in_len -= blk_len;
    blk_len = 5552;
  }
  return (unsigned)(s2 << 16) + (unsigned)s1;
}
extern int
zsdeflate(struct sdefl *s, void *out, const void *in, int n, int lvl) {
  int p = 0;
  unsigned a = 0;
  unsigned char *q = (unsigned char*)out;

  s->bits = s->bitcnt = 0;
  sdefl_put(&q, s, 0x78, 8); /* deflate, 32k window */
  sdefl_put(&q, s, 0x01, 8); /* fast compression */
  q += sdefl_compr(s, q, (const unsigned char*)in, n, lvl);

  /* append adler checksum */
  a = sdefl_adler32(SDEFL_ADLER_INIT, (const unsigned char*)in, n);
  for (p = 0; p < 4; ++p) {
    sdefl_put(&q, s, (a >> 24) & 0xFF, 8);
    a <<= 8;
  }
  return (int)(q - (unsigned char*)out);
}
extern int
sdefl_bound(int len) {
  int a = 128 + (len * 110) / 100;
  int b = 128 + len + ((len / (31 * 1024)) + 1) * 5;
  return (a > b) ? a : b;
}
#endif /* SDEFL_IMPLEMENTATION */

/*
# Small Deflate
`sdefl` is a small bare bone lossless compression library in ANSI C (ISO C90)
which implements the Deflate (RFC 1951) compressed data format specification standard.
It is mainly tuned to get as much speed and compression ratio from as little code
as needed to keep the implementation as concise as possible.

## Features
- Portable single header and source file duo written in ANSI C (ISO C90)
- Dual license with either MIT or public domain
- Small implementation
    - Deflate: 525 LoC
    - Inflate: 320 LoC
- Webassembly:
    - Deflate ~3.7 KB (~2.2KB compressed)
    - Inflate ~3.6 KB (~2.2KB compressed)

## Usage:
This file behaves differently depending on what symbols you define
before including it.

Header-File mode:
If you do not define `SINFL_IMPLEMENTATION` before including this file, it
will operate in header only mode. In this mode it declares all used structs
and the API of the library without including the implementation of the library.

Implementation mode:
If you define `SINFL_IMPLEMENTATION` before including this file, it will
compile the implementation. Make sure that you only include
this file implementation in *one* C or C++ file to prevent collisions.

### Benchmark

| Compressor name         | Compression| Decompress.| Compr. size | Ratio |
| ------------------------| -----------| -----------| ----------- | ----- |
| sdefl 1.0 -0            |   127 MB/s |   233 MB/s |    40004116 | 39.88 |
| sdefl 1.0 -1            |   111 MB/s |   259 MB/s |    38940674 | 38.82 |
| sdefl 1.0 -5            |    45 MB/s |   275 MB/s |    36577183 | 36.46 |
| sdefl 1.0 -7            |    38 MB/s |   276 MB/s |    36523781 | 36.41 |
| zlib 1.2.11 -1          |    72 MB/s |   307 MB/s |    42298774 | 42.30 |
| zlib 1.2.11 -6          |    24 MB/s |   313 MB/s |    36548921 | 36.55 |
| zlib 1.2.11 -9          |    20 MB/s |   314 MB/s |    36475792 | 36.48 |
| miniz 1.0 -1            |   122 MB/s |   208 MB/s |    48510028 | 48.51 |
| miniz 1.0 -6            |    27 MB/s |   260 MB/s |    36513697 | 36.51 |
| miniz 1.0 -9            |    23 MB/s |   261 MB/s |    36460101 | 36.46 |
| libdeflate 1.3 -1       |   147 MB/s |   667 MB/s |    39597378 | 39.60 |
| libdeflate 1.3 -6       |    69 MB/s |   689 MB/s |    36648318 | 36.65 |
| libdeflate 1.3 -9       |    13 MB/s |   672 MB/s |    35197141 | 35.20 |
| libdeflate 1.3 -12      |  8.13 MB/s |   670 MB/s |    35100568 | 35.10 |

### Compression
Results on the [Silesia compression corpus](http://sun.aei.polsl.pl/~sdeor/index.php?page=silesia):

| File    |   Original | `sdefl 0`      | `sdefl 5`     | `sdefl 7` |
| :------ | ---------: | -----------------: | ---------: | ----------: |
| dickens | 10.192.446 |  4,260,187|  3,845,261|   3,833,657 |
| mozilla | 51.220.480 | 20,774,706 | 19,607,009 |  19,565,867 |
| mr      |  9.970.564 | 3,860,531 |  3,673,460 |   3,665,627 |
| nci     | 33.553.445 | 4,030,283 |  3,094,526 |   3,006,075 |
| ooffice |  6.152.192 | 3,320,063 |  3,186,373 |   3,183,815 |
| osdb    | 10.085.684 | 3,919,646 |  3,649,510 |   3,649,477 |
| reymont |  6.627.202 | 2,263,378 |  1,857,588 |   1,827,237 |
| samba   | 21.606.400 | 6,121,797 |  5,462,670 |   5,450,762 |
| sao     |  7.251.944 | 5,612,421 |  5,485,380 |   5,481,765 |
| webster | 41.458.703 | 13,972,648 | 12,059,432 |  11,991,421 |
| xml     |  5.345.280 | 886,620|    674,009 |     662,141 |
| x-ray   |  8.474.240 | 6,304,655 |  6,244,779 |   6,244,779 |

## License
```
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2020 Micha Mettke
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
```
*/
#ifndef SINFL_H_INCLUDED
#define SINFL_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#define SINFL_PRE_TBL_SIZE 128
#define SINFL_LIT_TBL_SIZE 1334
#define SINFL_OFF_TBL_SIZE 402

struct sinfl {
  int bits, bitcnt;
  unsigned lits[SINFL_LIT_TBL_SIZE];
  unsigned dsts[SINFL_OFF_TBL_SIZE];
};
extern int sinflate(void *out, const void *in, int size);
extern int zsinflate(void *out, const void *in, int size);

#ifdef __cplusplus
}
#endif

#endif /* SINFL_H_INCLUDED */

#ifdef SINFL_IMPLEMENTATION

#include <string.h> /* memcpy, memset */

static int
sinfl_bsr(unsigned n) {
#ifdef _MSC_VER
  _BitScanReverse(&n, n);
  return n;
#elif defined(__GNUC__) || defined(__clang__)
  return 31 - __builtin_clz(n);
#endif
}
static int
sinfl_get(const unsigned char **src, const unsigned char *end, struct sinfl *s,
          int n) {
  const unsigned char *in = *src;
  int v = s->bits & ((1 << n)-1);
  s->bits >>= n;
  s->bitcnt = s->bitcnt - n;
  s->bitcnt = s->bitcnt < 0 ? 0 : s->bitcnt;
  while (s->bitcnt < 16 && in < end) {
    s->bits |= (*in++) << s->bitcnt;
    s->bitcnt += 8;
  }
  *src = in;
  return v;
}
struct sinfl_gen {
  int len;
  int cnt;
  int word;
  short* sorted;
};
static int
sinfl_build_tbl(struct sinfl_gen *gen, unsigned *tbl, int tbl_bits,
                const int *cnt) {
  int tbl_end = 0;
  while (!(gen->cnt = cnt[gen->len])) {
    ++gen->len;
  }
  tbl_end = 1 << gen->len;
  while (gen->len <= tbl_bits) {
    do {unsigned bit = 0;
      tbl[gen->word] = (*gen->sorted++ << 16) | gen->len;
      if (gen->word == tbl_end - 1) {
        for (; gen->len < tbl_bits; gen->len++) {
          memcpy(&tbl[tbl_end], tbl, (size_t)tbl_end * sizeof(tbl[0]));
          tbl_end <<= 1;
        }
        return 1;
      }
      bit = 1 << sinfl_bsr((unsigned)(gen->word ^ (tbl_end - 1)));
      gen->word &= bit - 1;
      gen->word |= bit;
    } while (--gen->cnt);
    do {
      if (++gen->len <= tbl_bits) {
        memcpy(&tbl[tbl_end], tbl, (size_t)tbl_end * sizeof(tbl[0]));
        tbl_end <<= 1;
      }
    } while (!(gen->cnt = cnt[gen->len]));
  }
  return 0;
}
static void
sinfl_build_subtbl(struct sinfl_gen *gen, unsigned *tbl, int tbl_bits,
                   const int *cnt) {
  int sub_bits = 0;
  int sub_start = 0;
  int sub_prefix = -1;
  int tbl_end = 1 << tbl_bits;
  while (1) {
    unsigned entry;
    int bit, stride, i;
    /* start new subtable */
    if ((gen->word & ((1 << tbl_bits)-1)) != sub_prefix) {
      int used = 0;
      sub_prefix = gen->word & ((1 << tbl_bits)-1);
      sub_start = tbl_end;
      sub_bits = gen->len - tbl_bits;
      used = gen->cnt;
      while (used < (1 << sub_bits)) {
        sub_bits++;
        used = (used << 1) + cnt[tbl_bits + sub_bits];
      }
      tbl_end = sub_start + (1 << sub_bits);
      tbl[sub_prefix] = (sub_start << 16) | 0x10 | (sub_bits & 0xf);
    }
    /* fill subtable */
    entry = (*gen->sorted << 16) | ((gen->len - tbl_bits) & 0xf);
    gen->sorted++;
    i = sub_start + (gen->word >> tbl_bits);
    stride = 1 << (gen->len - tbl_bits);
    do {
      tbl[i] = entry;
      i += stride;
    } while (i < tbl_end);
    if (gen->word == (1 << gen->len)-1) {
      return;
    }
    bit = 1 << sinfl_bsr(gen->word ^ ((1 << gen->len) - 1));
    gen->word &= bit - 1;
    gen->word |= bit;
    gen->cnt--;
    while (!gen->cnt) {
      gen->cnt = cnt[++gen->len];
    }
  }
}
static void
sinfl_build(unsigned *tbl, unsigned char *lens, int tbl_bits, int maxlen,
            int symcnt) {
  int i, used = 0;
  short sort[288];
  int cnt[16] = {0}, off[16]= {0};
  struct sinfl_gen gen = {0};
  gen.sorted = sort;
  gen.len = 1;

  for (i = 0; i < symcnt; ++i)
    cnt[lens[i]]++;
  off[1] = cnt[0];
  for (i = 1; i < maxlen; ++i) {
    off[i + 1] = off[i] + cnt[i];
    used = (used << 1) + cnt[i];
  }
  used = (used << 1) + cnt[i];
  for (i = 0; i < symcnt; ++i)
    gen.sorted[off[lens[i]]++] = (short)i;
  gen.sorted += off[0];

  if (used < (1 << maxlen)){
    for (i = 0; i < 1 << tbl_bits; ++i)
      tbl[i] = (0 << 16u) | 1;
    return;
  }
  if (!sinfl_build_tbl(&gen, tbl, tbl_bits, cnt)){
    sinfl_build_subtbl(&gen, tbl, tbl_bits, cnt);
  }
}
static int
sinfl_decode(const unsigned char **in, const unsigned char *end,
             struct sinfl *s, const unsigned *tbl, int bit_len) {
  int idx = s->bits & ((1 << bit_len) - 1);
  unsigned key = tbl[idx];
  if (key & 0x10) {
    /* sub-table lookup */
    int len = key & 0x0f;
    sinfl_get(in, end, s, bit_len);
    idx = s->bits & ((1 << len)-1);
    key = tbl[((key >> 16) & 0xffff) + (unsigned)idx];
  }
  sinfl_get(in, end, s, key & 0x0f);
  return (key >> 16) & 0x0fff;
}
static int
sinfl_decompress(unsigned char *out, const unsigned char *in, int size) {
  static const unsigned char order[] = {16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15};
  static const short dbase[30+2] = {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,
      257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};
  static const unsigned char dbits[30+2] = {0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,
      10,10,11,11,12,12,13,13,0,0};
  static const short lbase[29+2] = {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,
      43,51,59,67,83,99,115,131,163,195,227,258,0,0};
  static const unsigned char lbits[29+2] = {0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,
      4,4,4,5,5,5,5,0,0,0};

  const unsigned char *e = in + size, *o = out;
  enum sinfl_states {hdr,stored,fixed,dyn,blk};
  enum sinfl_states state = hdr;
  struct sinfl s = {0};
  int last = 0;

  sinfl_get(&in,e,&s,0); /* buffer input */
  while (in < e || s.bitcnt) {
    switch (state) {
    case hdr: {
      int type = 0; /* block header */
      last = sinfl_get(&in,e,&s,1);
      type = sinfl_get(&in,e,&s,2);

      switch (type) {default: return (int)(out-o);
      case 0x00: state = stored; break;
      case 0x01: state = fixed; break;
      case 0x02: state = dyn; break;}
    } break;
    case stored: {
      int len, nlen; /* uncompressed block */
      sinfl_get(&in,e,&s,s.bitcnt & 7);
      len = sinfl_get(&in,e,&s,16);
      nlen = sinfl_get(&in,e,&s,16);
      in -= 2; s.bitcnt = 0;

      if (len > (e-in) || !len)
        return (int)(out-o);
      memcpy(out, in, (size_t)len);
      in += len, out += len;
      state = hdr;
    } break;
    case fixed: {
      /* fixed huffman codes */
      int n; unsigned char lens[288+32];
      for (n = 0; n <= 143; n++) lens[n] = 8;
      for (n = 144; n <= 255; n++) lens[n] = 9;
      for (n = 256; n <= 279; n++) lens[n] = 7;
      for (n = 280; n <= 287; n++) lens[n] = 8;
      for (n = 0; n < 32; n++) lens[288+n] = 5;

      /* build lit/dist tables */
      sinfl_build(s.lits, lens, 10, 15, 288);
      sinfl_build(s.dsts, lens + 288, 8, 15, 32);
      state = blk;
    } break;
    case dyn: {
        /* dynamic huffman codes */
        int n, i;
        unsigned hlens[SINFL_PRE_TBL_SIZE];
        unsigned char nlens[19] = {0}, lens[288+32];
        int nlit = 257 + sinfl_get(&in,e,&s,5);
        int ndist = 1 + sinfl_get(&in,e,&s,5);
        int nlen = 4 + sinfl_get(&in,e,&s,4);
        for (n = 0; n < nlen; n++)
          nlens[order[n]] = (unsigned char)sinfl_get(&in,e,&s,3);
        sinfl_build(hlens, nlens, 7, 7, 19);

        /* decode code lengths */
        for (n = 0; n < nlit + ndist;) {
          int sym = sinfl_decode(&in, e, &s, hlens, 7);
          switch (sym) {default: lens[n++] = (unsigned char)sym; break;
          case 16: for (i=3+sinfl_get(&in,e,&s,2);i;i--,n++) lens[n]=lens[n-1]; break;
          case 17: for (i=3+sinfl_get(&in,e,&s,3);i;i--,n++) lens[n]=0; break;
          case 18: for (i=11+sinfl_get(&in,e,&s,7);i;i--,n++) lens[n]=0; break;}
        }
        /* build lit/dist tables */
        sinfl_build(s.lits, lens, 10, 15, nlit);
        sinfl_build(s.dsts, lens + nlit, 8, 15, ndist);
        state = blk;
    } break;
    case blk: {
      /* decompress block */
      int i, sym = sinfl_decode(&in, e, &s, s.lits, 10);
      if (sym > 256) {sym -= 257; /* match symbol */
        {int len = sinfl_get(&in, e, &s, lbits[sym]) + lbase[sym];
        int dsym = sinfl_decode(&in, e, &s, s.dsts, 8);
        int offs = sinfl_get(&in, e, &s, dbits[dsym]) + dbase[dsym];
        if (offs > (int)(out-o)) {
          return (int)(out-o);
        } else if (offs == 1) {
          /* rle match copying */
          unsigned char c = *(out - offs);
          unsigned long w = (c << 24) | (c << 16) | (c << 8) | c;
          for (i = 0; i < len >> 2; ++i) {
            memcpy(out, &w, 4);
            out += 4;
          }
          len = len & 3;
        } else if (offs >= 4) {
          /* copy match */
          int wcnt = len >> 2;
          for (i = 0; i < wcnt; ++i) {
            unsigned long w = 0;
            memcpy(&w, out - offs, 4);
            memcpy(out, &w, 4);
            out += 4;
          }
          len = len & 3;
        }
        for (i = 0; i < len; ++i)
          {*out = *(out-offs), out++;}
        }
      } else if (sym == 256) {
        /* end of block */
        if (last) return (int)(out-o);
        state = hdr;
        break;
        /* literal */
      } else *out++ = (unsigned char)sym;
    } break;}
  }
  return (int)(out-o);
}
extern int
sinflate(void *out, const void *in, int size) {
  return sinfl_decompress((unsigned char*)out, (const unsigned char*)in, size);
}
static unsigned
sinfl_adler32(unsigned adler32, const unsigned char *in, int in_len) {
  const unsigned ADLER_MOD = 65521;
  unsigned s1 = adler32 & 0xffff;
  unsigned s2 = adler32 >> 16;
  unsigned blk_len, i;

  blk_len = in_len % 5552;
  while (in_len) {
    for (i=0; i + 7 < blk_len; i += 8) {
      s1 += in[0]; s2 += s1;
      s1 += in[1]; s2 += s1;
      s1 += in[2]; s2 += s1;
      s1 += in[3]; s2 += s1;
      s1 += in[4]; s2 += s1;
      s1 += in[5]; s2 += s1;
      s1 += in[6]; s2 += s1;
      s1 += in[7]; s2 += s1;
      in += 8;
    }
    for (; i < blk_len; ++i)
      s1 += *in++, s2 += s1;
    s1 %= ADLER_MOD; s2 %= ADLER_MOD;
    in_len -= blk_len;
    blk_len = 5552;
  } return (unsigned)(s2 << 16) + (unsigned)s1;
}
extern int
zsinflate(void *out, const void *mem, int size) {
  const unsigned char *in = (const unsigned char*)mem;
  if (size >= 6) {
    const unsigned char *eob = in + size - 4;
    int n = sinfl_decompress((unsigned char*)out, in + 2u, size);
    unsigned a = sinfl_adler32(1u, (unsigned char*)out, n);
    unsigned h = eob[0] << 24 | eob[1] << 16 | eob[2] << 8 | eob[3] << 0;
    return a == h ? n : -1;
  } else {
    return -1;
  }
}
#endif

// end of deflate.c

unsigned deflate_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags) {
    struct sdefl *s = (struct sdefl *)malloc(sizeof(struct sdefl)); //{0};
    int rc = zsdeflate(s, (unsigned char *)out, (const unsigned char *)in, (int)inlen, flags > SDEFL_LVL_MAX ? SDEFL_LVL_MAX : flags);
    free(s);
    return rc > 0 ? (unsigned)rc : 0;
}
unsigned deflate_decode(const void *in, unsigned inlen, void *out, unsigned outlen) {
    int rc = sinflate((unsigned char *)out, (const unsigned char *)in + 2, (int)inlen - 2);
    return rc > 0 ? (unsigned)rc : 0;
}
unsigned deflate_bounds(unsigned inlen, unsigned flags) {
    return (unsigned)sdefl_bound(inlen);
}

#endif // DEFLATE_C


#ifdef DEFLATE_DEMO
#pragma once
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level=1;
    char out[128];
    unsigned outlen = deflate_encode(longcopy, strlen(longcopy)+1, out, 128, level);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    unsigned unpacked = deflate_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // DEFLATE_DEMO

#line 1 "amalgamated_lz4x.c" 
// LZ4X - An optimized LZ4 compressor
// Written and placed in the public domain by Ilya Muravyov (UNLICENSED)
// MemBased by @r-lyeh. @todo: thread-safe

unsigned lz4x_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags); //[1(fastest)..(6)..15(uber)]
unsigned lz4x_decode(const void *in, unsigned inlen, void *out, unsigned outlen);
unsigned lz4x_bounds(unsigned inlen, unsigned flags);


#ifdef LZ4X_C
#pragma once

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_DISABLE_PERFCRIT_LOCKS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

#ifndef LZ4X_REALLOC
#define LZ4X_REALLOC realloc
#endif

#define LZ4X_BLOCK_SIZE (8<<20) // 8 MB
#define LZ4X_PADDING_LITERALS 5

#define LZ4X_WINDOW_BITS 16
#define LZ4X_WINDOW_SIZE (1<<LZ4X_WINDOW_BITS)
#define LZ4X_WINDOW_MASK (LZ4X_WINDOW_SIZE-1)

#define LZ4X_MIN_MATCH 4

#define LZ4X_EXCESS (16+(LZ4X_BLOCK_SIZE/255))

#define LZ4X_MIN(a, b) (((a)<(b))?(a):(b))
#define LZ4X_MAX(a, b) (((a)>(b))?(a):(b))

#define LZ4X_LOAD_16(p) (*(const uint16_t*)(p))
#define LZ4X_LOAD_32(p) (*(const uint32_t*)(p))
#define LZ4X_STORE_16(p, x) (*(uint16_t*)(p)=(x))
#define LZ4X_COPY_32(d, s) (*(uint32_t*)(d)=LZ4X_LOAD_32(s))

#define LZ4X_HASH_BITS 18
#define LZ4X_HASH_SIZE (1<<LZ4X_HASH_BITS)
#define LZ4X_NIL (-1)

#define LZ4X_HASH_32(p) ((LZ4X_LOAD_32(p)*0x9E3779B9)>>(32-LZ4X_HASH_BITS))

#define lz4x_wild_copy(d,s,n) do { \
    LZ4X_COPY_32(d, s); \
    LZ4X_COPY_32(d+4, s+4); \
    for (int i_=8; i_<n; i_+=8) { \
        LZ4X_COPY_32(d+i_, s+i_); \
        LZ4X_COPY_32(d+4+i_, s+4+i_); \
    } \
} while(0)


int lz4x_compress_optimal(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen)
{
    static int head[LZ4X_HASH_SIZE];
    static int nodes[LZ4X_WINDOW_SIZE][2];
    struct lz4x_path
    {
        int cum;
        int len;
        int dist;
    } *path; //path[LZ4X_BLOCK_SIZE+1];
    path = (struct lz4x_path*)LZ4X_REALLOC(0, sizeof(path[0]) * (inlen+1));

    int n = (int)inlen;

    // Pass 1: Find all matches

    for (int i=0; i<LZ4X_HASH_SIZE; ++i)
        head[i]=LZ4X_NIL;

    for (int p=0; p<n; ++p)
    {
        int best_len=0;
        int dist=0;

        const int max_match=(n-LZ4X_PADDING_LITERALS)-p;
        if (max_match>=LZ4X_MAX(12-LZ4X_PADDING_LITERALS, LZ4X_MIN_MATCH))
        {
            const int limit=LZ4X_MAX(p-LZ4X_WINDOW_SIZE, LZ4X_NIL);

            int* left=&nodes[p&LZ4X_WINDOW_MASK][1];
            int* right=&nodes[p&LZ4X_WINDOW_MASK][0];

            int left_len=0;
            int right_len=0;

            const uint32_t h=LZ4X_HASH_32(&in[p]);
            int s=head[h];
            head[h]=p;

            while (s>limit)
            {
                int len=LZ4X_MIN(left_len, right_len);

                if (in[s+len]==in[p+len])
                {
                    while (++len<max_match && in[s+len]==in[p+len]);

                    if (len>best_len)
                    {
                        best_len=len;
                        dist=p-s;

                        if (len==max_match || len>=(1<<16))
                            break;
                    }
                }

                if (in[s+len]<in[p+len]) // in/out out/in ?
                {
                    *right=s;
                    right=&nodes[s&LZ4X_WINDOW_MASK][1];
                    s=*right;
                    right_len=len;
                }
                else
                {
                    *left=s;
                    left=&nodes[s&LZ4X_WINDOW_MASK][0];
                    s=*left;
                    left_len=len;
                }
            }

            *left=LZ4X_NIL;
            *right=LZ4X_NIL;
        }

        path[p].len=best_len;
        path[p].dist=dist;
    }

    // Pass 2: Build the shortest path

    path[n].cum=0;

    int count=15;

    for (int p=n-1; p>0; --p)
    {
        int c0=path[p+1].cum+1;

        if (--count==0)
        {
            count=255;
            ++c0;
        }

        int len=path[p].len;
        if (len>=LZ4X_MIN_MATCH)
        {
            int c1=1<<30;

            const int j=LZ4X_MAX(len-255, LZ4X_MIN_MATCH);
            for (int i=len; i>=j; --i)
            {
                int tmp=path[p+i].cum+3;

                if (i>=(15+LZ4X_MIN_MATCH))
                    tmp+=1+((i-(15+LZ4X_MIN_MATCH))/255);

                if (tmp<c1)
                {
                    c1=tmp;
                    len=i;
                }
            }

            if (c1<=c0)
            {
                path[p].cum=c1;
                path[p].len=len;

                count=15;
            }
            else
            {
                path[p].cum=c0;
                path[p].len=0;
            }
        }
        else
            path[p].cum=c0;
    }

    // Pass 3: Output the codes

    int op=0;
    int pp=0;

    int p=0;
    while (p<n)
    {
        if (path[p].len>=LZ4X_MIN_MATCH)
        {
            int len=path[p].len-LZ4X_MIN_MATCH;
            const int nib=LZ4X_MIN(len, 15);

            if (pp!=p)
            {
                const int run=p-pp;
                if (run>=15)
                {
                    out[op++]=(15<<4)+nib;

                    int j=run-15;
                    for (; j>=255; j-=255)
                        out[op++]=255;
                    out[op++]=j;
                }
                else
                    out[op++]=(run<<4)+nib;

                lz4x_wild_copy(&out[op], &in[pp], run);
                op+=run;
            }
            else
                out[op++]=nib;

            LZ4X_STORE_16(&out[op], path[p].dist);
            op+=2;

            if (len>=15)
            {
                len-=15;
                for (; len>=255; len-=255)
                    out[op++]=255;
                out[op++]=len;
            }

            p+=path[p].len;

            pp=p;
        }
        else
            ++p;
    }

    if (pp!=p)
    {
        const int run=p-pp;
        if (run>=15)
        {
            out[op++]=15<<4;

            int j=run-15;
            for (; j>=255; j-=255)
                out[op++]=255;
            out[op++]=j;
        }
        else
            out[op++]=run<<4;

        lz4x_wild_copy(&out[op], &in[pp], run);
        op+=run;
    }

    LZ4X_REALLOC(path, 0);

    const int comp_len=op;
    return comp_len;
}

int lz4x_compress(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen, unsigned max_chain)
{
    static int head[LZ4X_HASH_SIZE];
    static int tail[LZ4X_WINDOW_SIZE];

    int n = (int)inlen;

    for (int i=0; i<LZ4X_HASH_SIZE; ++i)
        head[i]=LZ4X_NIL;

    int op=0;
    int pp=0;

    int p=0;
    while (p<n)
    {
        int best_len=0;
        int dist=0;

        const int max_match=(n-LZ4X_PADDING_LITERALS)-p;
        if (max_match>=LZ4X_MAX(12-LZ4X_PADDING_LITERALS, LZ4X_MIN_MATCH))
        {
            const int limit=LZ4X_MAX(p-LZ4X_WINDOW_SIZE, LZ4X_NIL);
            int chain_len=max_chain;

            int s=head[LZ4X_HASH_32(&in[p])];
            while (s>limit)
            {
                if (in[s+best_len]==in[p+best_len] && LZ4X_LOAD_32(&in[s])==LZ4X_LOAD_32(&in[p]))
                {
                    int len=LZ4X_MIN_MATCH;
                    while (len<max_match && in[s+len]==in[p+len])
                        ++len;

                    if (len>best_len)
                    {
                        best_len=len;
                        dist=p-s;

                        if (len==max_match)
                            break;
                    }
                }

                if (--chain_len<=0)
                    break;

                s=tail[s&LZ4X_WINDOW_MASK];
            }
        }

        if (best_len>=LZ4X_MIN_MATCH)
        {
            int len=best_len-LZ4X_MIN_MATCH;
            const int nib=LZ4X_MIN(len, 15);

            if (pp!=p)
            {
                const int run=p-pp;
                if (run>=15)
                {
                    out[op++]=(15<<4)+nib;

                    int j=run-15;
                    for (; j>=255; j-=255)
                        out[op++]=255;
                    out[op++]=j;
                }
                else
                    out[op++]=(run<<4)+nib;

                lz4x_wild_copy(&out[op], &in[pp], run);
                op+=run;
            }
            else
                out[op++]=nib;

            LZ4X_STORE_16(&out[op], dist);
            op+=2;

            if (len>=15)
            {
                len-=15;
                for (; len>=255; len-=255)
                    out[op++]=255;
                out[op++]=len;
            }

            pp=p+best_len;

            while (p<pp)
            {
                const uint32_t h=LZ4X_HASH_32(&in[p]); // out?
                tail[p&LZ4X_WINDOW_MASK]=head[h];
                head[h]=p++;
            }
        }
        else
        {
            const uint32_t h=LZ4X_HASH_32(&in[p]); // out?
            tail[p&LZ4X_WINDOW_MASK]=head[h];
            head[h]=p++;
        }
    }

    if (pp!=p)
    {
        const int run=p-pp;
        if (run>=15)
        {
            out[op++]=15<<4;

            int j=run-15;
            for (; j>=255; j-=255)
                out[op++]=255;
            out[op++]=j;
        }
        else
            out[op++]=run<<4;

        lz4x_wild_copy(&out[op], &in[pp], run);
        op+=run;
    }

    const int comp_len=op;
    return comp_len;
}

int lz4x_decompress(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen)
{
    int n = (int)inlen;

    int p=0;
    int ip=0;
    const int ip_end=ip+n;

    for (;;)
    {
        const int token=in[ip++];
        if (token>=16)
        {
            int run=token>>4;
            if (run==15)
            {
                for (;;)
                {
                    const int c=in[ip++];
                    run+=c;
                    if (c!=255)
                        break;
                }
            }
            if ((p+run)>outlen) return 0; // -1

            lz4x_wild_copy(&out[p], &in[ip], run);
            p+=run;
            ip+=run;
            if (ip>=ip_end)
                break;
        }

        int s=p-LZ4X_LOAD_16(&in[ip]);
        ip+=2;
        if (s<0) return 0; // -1

        int len=(token&15)+LZ4X_MIN_MATCH;
        if (len==(15+LZ4X_MIN_MATCH))
        {
            for (;;)
            {
                const int c=in[ip++];
                len+=c;
                if (c!=255)
                    break;
            }
        }
        if ((p+len)>outlen) return 0; // -1

        if ((p-s)>=4)
        {
            lz4x_wild_copy(&out[p], &out[s], len);
            p+=len;
        }
        else
        {
            while (len--!=0)
                out[p++]=out[s++];
        }
    }

    return p;
}

unsigned lz4x_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags/*[1..15*/) {
    unsigned level = (unsigned)(flags > 15 ? 15 : flags < 1 ? 1 : flags);
    if(level >= 15) return lz4x_compress_optimal((const uint8_t*)in, inlen, (uint8_t*)out, outlen);
    return (unsigned)lz4x_compress((const uint8_t*)in, inlen, (uint8_t*)out, outlen, level);
}
unsigned lz4x_decode(const void *in, unsigned inlen, void *out, unsigned outlen) {
    return (unsigned)lz4x_decompress((const uint8_t*)in, (size_t)inlen, (uint8_t*)out, (size_t)outlen);
}
unsigned lz4x_bounds(unsigned inlen, unsigned flags) {
    return (unsigned)(inlen + (inlen/255) + 16);
}

#endif // LZ4X_C


#ifdef LZ4X_DEMO
#pragma once
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level=1;
    char out[128];
    unsigned outlen = lz4x_encode(longcopy, strlen(longcopy)+1, out, 128, level);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    unsigned unpacked = lz4x_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // LZ4X_DEMO

#line 1 "amalgamated_lzma.c" 
// LzFind.c  -- Match finder for LZ algorithms 2009-04-22 : Igor Pavlov : Public domain
// LzmaDec.c -- LZMA Decoder                   2009-09-20 : Igor Pavlov : Public domain
// LzmaEnc.c -- LZMA Encoder                   2009-11-24 : Igor Pavlov : Public domain
// Additional code by @r-lyeh, public domain. TOC: glue.h+lzfind.h/c+lzmaenc.h/c+lzmadec.h/c+glue.c

unsigned lzma_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags); // [0..(7)..9]
unsigned lzma_decode(const void *in, unsigned inlen, void *out, unsigned outlen);
unsigned lzma_bounds(unsigned inlen, unsigned flags);


#ifdef LZMA_C
#pragma once

// glue.h

#ifndef LZMA_REALLOC
#define LZMA_REALLOC realloc
#endif

#define LZMA_MALLOC(s) LZMA_REALLOC(0, s)
#define LZMA_FREE(p)   LZMA_REALLOC(p, 0)

#define _FILE_OFFSET_BITS 64
#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#ifndef max
#define max(x,y) ((x) >= (y) ? (x) : (y))
#endif
#ifndef min
#define min(x,y) ((x) <= (y) ? (x) : (y))
#endif

#ifdef _WIN32
//#include <io.h>
#else
//#include <unistd.h>
#endif


/* #define SHOW_STAT */
/* #define SHOW_STAT2 */

typedef int State;

enum {
        min_dictionary_bits = 12,
        min_dictionary_size = 1 << min_dictionary_bits,
        max_dictionary_bits = 29,
        max_dictionary_size = 1 << max_dictionary_bits,
        max_dictionary_bits_c = 27,           /* kDicLogSizeMaxCompress */
        max_dictionary_size_c = 1 << max_dictionary_bits_c,
        literal_context_bits = 3,
        literal_pos_state_bits = 0,               /* not used */
        pos_state_bits = 2,

        len_low_bits = 3,
        len_mid_bits = 3,
        len_high_bits = 8,
        len_low_symbols = 1 << len_low_bits,
        len_mid_symbols = 1 << len_mid_bits,
        len_high_symbols = 1 << len_high_bits,
        max_len_symbols = len_low_symbols + len_mid_symbols + len_high_symbols,

        min_match_len = 2,                    /* must be 2 */
        max_match_len = min_match_len + max_len_symbols - 1,  /* 273 */
        min_match_len_limit = 5
};

enum { 
        SZ_OK = 0, 
        SZ_ERROR_READ = 8, 
        SZ_ERROR_WRITE = 9,
};

// io interface
static int readblock( const int fd, uint8_t *buf,int size );
static int writeblock( const int fd, const uint8_t *buf, int size );

/* LzFind.h -- Match finder for LZ algorithms
2009-04-22 : Igor Pavlov : Public domain */

typedef uint32_t CLzRef;

typedef struct
{
    uint8_t *bufferBase;
    uint8_t *buffer;
    CLzRef *hash;
    CLzRef *son;
    uint32_t pos;
    uint32_t posLimit;
    uint32_t streamPos;
    uint32_t lenLimit;

    uint32_t cyclicBufferPos;
    uint32_t cyclicBufferSize; /* it must be = (historySize + 1) */

    uint32_t matchMaxLen;
    uint32_t hashMask;
    uint32_t cutValue;

    uint32_t blockSize;
    uint32_t keepSizeBefore;
    uint32_t keepSizeAfter;

    uint32_t numHashBytes;
    uint32_t historySize;
    uint32_t hashSizeSum;
    uint32_t numSons;
    int infd;
    int result;
    uint32_t crc;
    bool btMode;
    bool streamEndWasReached;
} CMatchFinder;


/* Conditions:
         historySize <= 3 GB
         keepAddBufferBefore + matchMaxLen + keepAddBufferAfter < 511MB
*/
int Mf_Init(CMatchFinder *p, const int ifd, const int mc, uint32_t historySize,
        uint32_t keepAddBufferBefore, uint32_t matchMaxLen, uint32_t keepAddBufferAfter);

void Mf_Free(CMatchFinder *p);


/*
Conditions:
    Mf_GetNumAvailableBytes_Func must be called before each Mf_GetMatchLen_Func.
    Mf_GetPointerToCurrentPos_Func's result must be used only before any other function
*/

typedef uint32_t (*Mf_GetMatches_Func)(void *object, uint32_t *distances);
typedef void (*Mf_Skip_Func)(void *object, uint32_t);

typedef struct _IMatchFinder
{
    Mf_GetMatches_Func GetMatches;
    Mf_Skip_Func Skip;
} IMatchFinder;

void Mf_CreateVTable(CMatchFinder *p, IMatchFinder *vTable);

static inline uint32_t Mf_GetNumAvailableBytes(CMatchFinder *p)
    { return p->streamPos - p->pos; }

static inline uint8_t Mf_GetIndexByte(CMatchFinder *p, int index)
    { return p->buffer[index]; }

static inline uint8_t * Mf_GetPointerToCurrentPos(CMatchFinder *p)
    { return p->buffer; }

/* LzFind.c -- Match finder for LZ algorithms
2009-04-22 : Igor Pavlov : Public domain */

static uint32_t crc32[256];    /* Table of CRCs of all 8-bit messages. */

static inline void CRC32_init(void) {
        for( unsigned n = 0; n < 256; ++n ) {
                unsigned c = n;
                for( int k = 0; k < 8; ++k ) {
                        if( c & 1 ) c = 0xEDB88320U ^ ( c >> 1 ); else c >>= 1;
                }
                crc32[n] = c;
        }
}

static inline void CRC32_update_buf(uint32_t* const crc, const uint8_t* const buffer, const int size) {
        uint32_t c = *crc;
        for( int i = 0; i < size; ++i )
                c = crc32[(c^buffer[i])&0xFF] ^ ( c >> 8 );
        *crc = c;
}

#define kHash2Size (1 << 10)
#define kHash3Size (1 << 16)
#define kHash4Size (1 << 20)

#define kFix3HashSize (kHash2Size)
#define kFix4HashSize (kHash2Size + kHash3Size)

#define HASH2_CALC hashValue = cur[0] | ((uint32_t)cur[1] << 8);

#define HASH3_CALC { \
    uint32_t temp = crc32[cur[0]] ^ cur[1]; \
    hash2Value = temp & (kHash2Size - 1); \
    hashValue = (temp ^ ((uint32_t)cur[2] << 8)) & p->hashMask; }

#define HASH4_CALC { \
    uint32_t temp = crc32[cur[0]] ^ cur[1]; \
    hash2Value = temp & (kHash2Size - 1); \
    hash3Value = (temp ^ ((uint32_t)cur[2] << 8)) & (kHash3Size - 1); \
    hashValue = (temp ^ ((uint32_t)cur[2] << 8) ^ (crc32[cur[3]] << 5)) & p->hashMask; }

#define kEmptyHashValue 0
#define kMaxValForNormalize ((uint32_t)0xFFFFFFFF)
#define kNormalizeStepMin (1 << 10) /* it must be power of 2 */
#define kNormalizeMask (~(kNormalizeStepMin - 1))

#define kStartMaxLen 3


static void Mf_ReadBlock(CMatchFinder *p)
{
    if (p->streamEndWasReached || p->result != SZ_OK)
        return;
    for (;;)
    {
        uint8_t * const dest = p->buffer + (p->streamPos - p->pos);
        const int size = (p->bufferBase + p->blockSize - dest);
        int rd;
        if (size == 0)
            return;
        rd = readblock( p->infd, dest, size );
        if (rd != size && errno)
            { p->result = SZ_ERROR_READ; return; }
        if (rd == 0)
        {
            p->streamEndWasReached = true;
            return;
        }
        CRC32_update_buf( &p->crc, dest, rd );
        p->streamPos += rd;
        if (p->streamPos - p->pos > p->keepSizeAfter)
            return;
    }
}


static void Mf_CheckAndMoveAndRead(CMatchFinder *p)
{
    if ((uint32_t)(p->bufferBase + p->blockSize - p->buffer) <= p->keepSizeAfter)
        {
        memmove(p->bufferBase,
            p->buffer - p->keepSizeBefore,
            p->streamPos - p->pos + p->keepSizeBefore);
        p->buffer = p->bufferBase + p->keepSizeBefore;
        }
    Mf_ReadBlock(p);
}


void Mf_Free(CMatchFinder *p)
{
    LZMA_FREE(p->hash);
    p->hash = 0;
    LZMA_FREE(p->bufferBase);
    p->bufferBase = 0;
}

static CLzRef* AllocRefs(uint32_t num)
{
    uint32_t sizeInBytes = num * sizeof(CLzRef);
    if (sizeInBytes / sizeof(CLzRef) != num)
        return 0;
    return (CLzRef *)LZMA_MALLOC(sizeInBytes);
}

static void Mf_SetLimits(CMatchFinder *p)
{
    uint32_t limit = kMaxValForNormalize - p->pos;
    uint32_t limit2 = p->cyclicBufferSize - p->cyclicBufferPos;
    if (limit2 < limit)
        limit = limit2;
    limit2 = p->streamPos - p->pos;
    if (limit2 <= p->keepSizeAfter)
    {
        if (limit2 > 0)
            limit2 = 1;
    }
    else
        limit2 -= p->keepSizeAfter;
    if (limit2 < limit)
        limit = limit2;
    {
        uint32_t lenLimit = p->streamPos - p->pos;
        if (lenLimit > p->matchMaxLen)
            lenLimit = p->matchMaxLen;
        p->lenLimit = lenLimit;
    }
    p->posLimit = p->pos + limit;
}


int Mf_Init(CMatchFinder *p, const int ifd, const int mc, uint32_t historySize,
        uint32_t keepAddBufferBefore, uint32_t matchMaxLen, uint32_t keepAddBufferAfter)
{
    const uint32_t sizeReserv = ( historySize >> 1 ) +
        (keepAddBufferBefore + matchMaxLen + keepAddBufferAfter) / 2 + (1 << 19);

    p->hash = 0;
    p->cutValue = mc;
    p->infd = ifd;
    p->btMode = true;
    p->numHashBytes = 4;
    p->crc = 0xFFFFFFFFU;
    p->keepSizeBefore = historySize + keepAddBufferBefore + 1;
    p->keepSizeAfter = matchMaxLen + keepAddBufferAfter;
    /* we need one additional byte, since we use MoveBlock after pos++ and before dictionary using */
    /* keepSizeBefore + keepSizeAfter + sizeReserv must be < 4G) */
    p->blockSize = p->keepSizeBefore + p->keepSizeAfter + sizeReserv;
    p->buffer = p->bufferBase = (uint8_t *)LZMA_MALLOC(p->blockSize);
    if( p->bufferBase )
    {
        uint32_t newCyclicBufferSize = historySize + 1;
        uint32_t hs;
        p->matchMaxLen = matchMaxLen;
        {
            if (p->numHashBytes == 2)
                hs = (1 << 16) - 1;
            else
            {
                hs = historySize - 1;
                hs |= (hs >> 1);
                hs |= (hs >> 2);
                hs |= (hs >> 4);
                hs |= (hs >> 8);
                hs >>= 1;
                hs |= 0xFFFF; /* don't change it! It's required for Deflate */
                if (hs > (1 << 24))
                {
                    if (p->numHashBytes == 3)
                        hs = (1 << 24) - 1;
                    else
                        hs >>= 1;
                }
            }
            p->hashMask = hs;
            hs++;
            if (p->numHashBytes > 2) hs += kHash2Size;
            if (p->numHashBytes > 3) hs += kHash3Size;
            if (p->numHashBytes > 4) hs += kHash4Size;
        }

        {
            uint32_t newSize;
            p->historySize = historySize;
            p->hashSizeSum = hs;
            p->cyclicBufferSize = newCyclicBufferSize;
            p->numSons = (p->btMode ? newCyclicBufferSize * 2 : newCyclicBufferSize);
            newSize = p->hashSizeSum + p->numSons;
            p->hash = AllocRefs(newSize);
            if (p->hash != 0)
            {
                uint32_t i;
                p->son = p->hash + p->hashSizeSum;
                for (i = 0; i < p->hashSizeSum; i++)
                    p->hash[i] = kEmptyHashValue;
                p->cyclicBufferPos = 0;
                p->pos = p->streamPos = p->cyclicBufferSize;
                p->result = SZ_OK;
                p->streamEndWasReached = false;
                Mf_ReadBlock(p);
                Mf_SetLimits(p);
                return 1;
            }
        }
    }
    Mf_Free(p);
    return 0;
}

static void Mf_Normalize3(uint32_t subValue, CLzRef *items, uint32_t numItems)
{
    uint32_t i;
    for (i = 0; i < numItems; i++)
    {
        uint32_t value = items[i];
        if (value <= subValue)
            value = kEmptyHashValue;
        else
            value -= subValue;
        items[i] = value;
    }
}

static void Mf_Normalize(CMatchFinder *p)
{
    uint32_t subValue = (p->pos - p->historySize - 1) & kNormalizeMask;
    Mf_Normalize3(subValue, p->hash, p->hashSizeSum + p->numSons);
    p->posLimit -= subValue;
    p->pos -= subValue;
    p->streamPos -= subValue;
}

static void Mf_CheckLimits(CMatchFinder *p)
{
    if (p->pos == kMaxValForNormalize)
        Mf_Normalize(p);
    if (!p->streamEndWasReached && p->keepSizeAfter == p->streamPos - p->pos)
        Mf_CheckAndMoveAndRead(p);
    if (p->cyclicBufferPos == p->cyclicBufferSize)
        p->cyclicBufferPos = 0;
    Mf_SetLimits(p);
}

static uint32_t * Hc_GetMatchesSpec(uint32_t lenLimit, uint32_t curMatch, uint32_t pos, const uint8_t *cur, CLzRef *son,
        uint32_t _cyclicBufferPos, uint32_t _cyclicBufferSize, uint32_t cutValue,
        uint32_t *distances, uint32_t maxLen)
{
    son[_cyclicBufferPos] = curMatch;
    for (;;)
    {
        uint32_t delta = pos - curMatch;
        if (cutValue-- == 0 || delta >= _cyclicBufferSize)
            return distances;
        {
            const uint8_t *pb = cur - delta;
            curMatch = son[_cyclicBufferPos - delta + ((delta > _cyclicBufferPos) ? _cyclicBufferSize : 0)];
            if (pb[maxLen] == cur[maxLen] && *pb == *cur)
            {
                uint32_t len = 0;
                while (++len != lenLimit)
                    if (pb[len] != cur[len])
                        break;
                if (maxLen < len)
                {
                    *distances++ = maxLen = len;
                    *distances++ = delta - 1;
                    if (len == lenLimit)
                        return distances;
                }
            }
        }
    }
}


static uint32_t * GetMatchesSpec1( uint32_t lenLimit, uint32_t curMatch,
        uint32_t pos, const uint8_t *cur, CLzRef *son,
        uint32_t _cyclicBufferPos, uint32_t _cyclicBufferSize, uint32_t cutValue,
        uint32_t *distances, uint32_t maxLen )
{
    CLzRef *ptr0 = son + (_cyclicBufferPos << 1) + 1;
    CLzRef *ptr1 = son + (_cyclicBufferPos << 1);
    uint32_t len0 = 0, len1 = 0;
    for (;;)
    {
        uint32_t delta = pos - curMatch;
        if (cutValue-- == 0 || delta >= _cyclicBufferSize)
        {
            *ptr0 = *ptr1 = kEmptyHashValue;
            return distances;
        }
        {
            CLzRef *pair = son + ((_cyclicBufferPos - delta + ((delta > _cyclicBufferPos) ? _cyclicBufferSize : 0)) << 1);
            const uint8_t *pb = cur - delta;
            uint32_t len = (len0 < len1 ? len0 : len1);
            if (pb[len] == cur[len])
            {
                if (++len != lenLimit && pb[len] == cur[len])
                    while (++len != lenLimit)
                        if (pb[len] != cur[len])
                            break;
                if (maxLen < len)
                {
                    *distances++ = maxLen = len;
                    *distances++ = delta - 1;
                    if (len == lenLimit)
                    {
                        *ptr1 = pair[0];
                        *ptr0 = pair[1];
                        return distances;
                    }
                }
            }
            if (pb[len] < cur[len])
            {
                *ptr1 = curMatch;
                ptr1 = pair + 1;
                curMatch = *ptr1;
                len1 = len;
            }
            else
            {
                *ptr0 = curMatch;
                ptr0 = pair;
                curMatch = *ptr0;
                len0 = len;
            }
        }
    }
}

static void SkipMatchesSpec(uint32_t lenLimit, uint32_t curMatch, uint32_t pos, const uint8_t *cur, CLzRef *son,
        uint32_t _cyclicBufferPos, uint32_t _cyclicBufferSize, uint32_t cutValue)
{
    CLzRef *ptr0 = son + (_cyclicBufferPos << 1) + 1;
    CLzRef *ptr1 = son + (_cyclicBufferPos << 1);
    uint32_t len0 = 0, len1 = 0;
    for (;;)
    {
        uint32_t delta = pos - curMatch;
        if (cutValue-- == 0 || delta >= _cyclicBufferSize)
        {
            *ptr0 = *ptr1 = kEmptyHashValue;
            return;
        }
        {
            CLzRef *pair = son + ((_cyclicBufferPos - delta + ((delta > _cyclicBufferPos) ? _cyclicBufferSize : 0)) << 1);
            const uint8_t *pb = cur - delta;
            uint32_t len = (len0 < len1 ? len0 : len1);
            if (pb[len] == cur[len])
            {
                while (++len != lenLimit)
                    if (pb[len] != cur[len])
                        break;
                {
                    if (len == lenLimit)
                    {
                        *ptr1 = pair[0];
                        *ptr0 = pair[1];
                        return;
                    }
                }
            }
            if (pb[len] < cur[len])
            {
                *ptr1 = curMatch;
                ptr1 = pair + 1;
                curMatch = *ptr1;
                len1 = len;
            }
            else
            {
                *ptr0 = curMatch;
                ptr0 = pair;
                curMatch = *ptr0;
                len0 = len;
            }
        }
    }
}

#define MOVE_POS \
    ++p->cyclicBufferPos; \
    p->buffer++; \
    if (++p->pos == p->posLimit) Mf_CheckLimits(p);

#define MOVE_POS_RET MOVE_POS return offset;

static void Mf_MovePos(CMatchFinder *p) { MOVE_POS; }

#define GET_MATCHES_HEADER2(minLen, ret_op) \
    uint32_t lenLimit; uint32_t hashValue; const uint8_t *cur; uint32_t curMatch; \
    lenLimit = p->lenLimit; { if (lenLimit < minLen) { Mf_MovePos(p); ret_op; }} \
    cur = p->buffer;

#define GET_MATCHES_HEADER(minLen) GET_MATCHES_HEADER2(minLen, return 0)
#define SKIP_HEADER(minLen)        GET_MATCHES_HEADER2(minLen, continue)

#define MF_PARAMS(p) p->pos, p->buffer, p->son, p->cyclicBufferPos, p->cyclicBufferSize, p->cutValue

#define GET_MATCHES_FOOTER(offset, maxLen) \
    offset = (uint32_t)(GetMatchesSpec1(lenLimit, curMatch, MF_PARAMS(p), \
    distances + offset, maxLen) - distances); MOVE_POS_RET;

#define SKIP_FOOTER \
    SkipMatchesSpec(lenLimit, curMatch, MF_PARAMS(p)); MOVE_POS;

static uint32_t Bt2_MatchFinder_GetMatches(CMatchFinder *p, uint32_t *distances)
{
    uint32_t offset;
    GET_MATCHES_HEADER(2)
    HASH2_CALC;
    curMatch = p->hash[hashValue];
    p->hash[hashValue] = p->pos;
    offset = 0;
    GET_MATCHES_FOOTER(offset, 1)
}

static uint32_t Bt3_MatchFinder_GetMatches(CMatchFinder *p, uint32_t *distances)
{
    uint32_t hash2Value, delta2, maxLen, offset;
    GET_MATCHES_HEADER(3)

    HASH3_CALC;

    delta2 = p->pos - p->hash[hash2Value];
    curMatch = p->hash[kFix3HashSize + hashValue];

    p->hash[hash2Value] =
    p->hash[kFix3HashSize + hashValue] = p->pos;


    maxLen = 2;
    offset = 0;
    if (delta2 < p->cyclicBufferSize && *(cur - delta2) == *cur)
    {
        for (; maxLen != lenLimit; maxLen++)
            if (cur[(ptrdiff_t)maxLen - delta2] != cur[maxLen])
                break;
        distances[0] = maxLen;
        distances[1] = delta2 - 1;
        offset = 2;
        if (maxLen == lenLimit)
        {
            SkipMatchesSpec(lenLimit, curMatch, MF_PARAMS(p));
            MOVE_POS_RET;
        }
    }
    GET_MATCHES_FOOTER(offset, maxLen)
}

static uint32_t Bt4_MatchFinder_GetMatches(CMatchFinder *p, uint32_t *distances)
{
    uint32_t hash2Value, hash3Value, delta2, delta3, maxLen, offset;
    GET_MATCHES_HEADER(4)

    HASH4_CALC;

    delta2 = p->pos - p->hash[                hash2Value];
    delta3 = p->pos - p->hash[kFix3HashSize + hash3Value];
    curMatch = p->hash[kFix4HashSize + hashValue];

    p->hash[                hash2Value] =
    p->hash[kFix3HashSize + hash3Value] =
    p->hash[kFix4HashSize + hashValue] = p->pos;

    maxLen = 1;
    offset = 0;
    if (delta2 < p->cyclicBufferSize && *(cur - delta2) == *cur)
    {
        distances[0] = maxLen = 2;
        distances[1] = delta2 - 1;
        offset = 2;
    }
    if (delta2 != delta3 && delta3 < p->cyclicBufferSize && *(cur - delta3) == *cur)
    {
        maxLen = 3;
        distances[offset + 1] = delta3 - 1;
        offset += 2;
        delta2 = delta3;
    }
    if (offset != 0)
    {
        for (; maxLen != lenLimit; maxLen++)
            if (cur[(ptrdiff_t)maxLen - delta2] != cur[maxLen])
                break;
        distances[offset - 2] = maxLen;
        if (maxLen == lenLimit)
        {
            SkipMatchesSpec(lenLimit, curMatch, MF_PARAMS(p));
            MOVE_POS_RET;
        }
    }
    if (maxLen < 3)
        maxLen = 3;
    GET_MATCHES_FOOTER(offset, maxLen)
}

static uint32_t Hc4_MatchFinder_GetMatches(CMatchFinder *p, uint32_t *distances)
{
    uint32_t hash2Value, hash3Value, delta2, delta3, maxLen, offset;
    GET_MATCHES_HEADER(4)

    HASH4_CALC;

    delta2 = p->pos - p->hash[                hash2Value];
    delta3 = p->pos - p->hash[kFix3HashSize + hash3Value];
    curMatch = p->hash[kFix4HashSize + hashValue];

    p->hash[                hash2Value] =
    p->hash[kFix3HashSize + hash3Value] =
    p->hash[kFix4HashSize + hashValue] = p->pos;

    maxLen = 1;
    offset = 0;
    if (delta2 < p->cyclicBufferSize && *(cur - delta2) == *cur)
    {
        distances[0] = maxLen = 2;
        distances[1] = delta2 - 1;
        offset = 2;
    }
    if (delta2 != delta3 && delta3 < p->cyclicBufferSize && *(cur - delta3) == *cur)
    {
        maxLen = 3;
        distances[offset + 1] = delta3 - 1;
        offset += 2;
        delta2 = delta3;
    }
    if (offset != 0)
    {
        for (; maxLen != lenLimit; maxLen++)
            if (cur[(ptrdiff_t)maxLen - delta2] != cur[maxLen])
                break;
        distances[offset - 2] = maxLen;
        if (maxLen == lenLimit)
        {
            p->son[p->cyclicBufferPos] = curMatch;
            MOVE_POS_RET;
        }
    }
    if (maxLen < 3)
        maxLen = 3;
    offset = (uint32_t)(Hc_GetMatchesSpec(lenLimit, curMatch, MF_PARAMS(p),
        distances + offset, maxLen) - (distances));
    MOVE_POS_RET
}


static void Bt2_MatchFinder_Skip(CMatchFinder *p, uint32_t num)
{
    do
    {
        SKIP_HEADER(2)
        HASH2_CALC;
        curMatch = p->hash[hashValue];
        p->hash[hashValue] = p->pos;
        SKIP_FOOTER
    }
    while (--num != 0);
}


static void Bt3_MatchFinder_Skip(CMatchFinder *p, uint32_t num)
{
    do
    {
        uint32_t hash2Value;
        SKIP_HEADER(3)
        HASH3_CALC;
        curMatch = p->hash[kFix3HashSize + hashValue];
        p->hash[hash2Value] =
        p->hash[kFix3HashSize + hashValue] = p->pos;
        SKIP_FOOTER
    }
    while (--num != 0);
}

static void Bt4_MatchFinder_Skip(CMatchFinder *p, uint32_t num)
{
    do
    {
        uint32_t hash2Value, hash3Value;
        SKIP_HEADER(4)
        HASH4_CALC;
        curMatch = p->hash[kFix4HashSize + hashValue];
        p->hash[                hash2Value] =
        p->hash[kFix3HashSize + hash3Value] = p->pos;
        p->hash[kFix4HashSize + hashValue] = p->pos;
        SKIP_FOOTER
    }
    while (--num != 0);
}

static void Hc4_MatchFinder_Skip(CMatchFinder *p, uint32_t num)
{
    do
    {
        uint32_t hash2Value, hash3Value;
        SKIP_HEADER(4)
        HASH4_CALC;
        curMatch = p->hash[kFix4HashSize + hashValue];
        p->hash[                hash2Value] =
        p->hash[kFix3HashSize + hash3Value] =
        p->hash[kFix4HashSize + hashValue] = p->pos;
        p->son[p->cyclicBufferPos] = curMatch;
        MOVE_POS
    }
    while (--num != 0);
}


void Mf_CreateVTable(CMatchFinder *p, IMatchFinder *vTable)
{
    if (!p->btMode)
    {
        vTable->GetMatches = (Mf_GetMatches_Func)Hc4_MatchFinder_GetMatches;
        vTable->Skip = (Mf_Skip_Func)Hc4_MatchFinder_Skip;
    }
    else if (p->numHashBytes == 2)
    {
        vTable->GetMatches = (Mf_GetMatches_Func)Bt2_MatchFinder_GetMatches;
        vTable->Skip = (Mf_Skip_Func)Bt2_MatchFinder_Skip;
    }
    else if (p->numHashBytes == 3)
    {
        vTable->GetMatches = (Mf_GetMatches_Func)Bt3_MatchFinder_GetMatches;
        vTable->Skip = (Mf_Skip_Func)Bt3_MatchFinder_Skip;
    }
    else
    {
        vTable->GetMatches = (Mf_GetMatches_Func)Bt4_MatchFinder_GetMatches;
        vTable->Skip = (Mf_Skip_Func)Bt4_MatchFinder_Skip;
    }
}

/* LzmaEnc.h -- LZMA Encoder
2009-02-07 : Igor Pavlov : Public domain */

/* ---------- CLzmaEncHandle Interface ---------- */

/* LzmaEnc_* functions can return the following exit codes:
Returns:
    SZ_OK           - OK
    SZ_ERROR_WRITE  - Write callback error.
*/

typedef void * CLzmaEncHandle;

CLzmaEncHandle LzmaEnc_Init( const int dict_size, const int match_len_limit,
                                                         const int infd, const int outfd );
void LzmaEnc_Free(CLzmaEncHandle p);
int LzmaEnc_Encode(CLzmaEncHandle p);

/* LzmaEnc.c -- LZMA Encoder
2009-11-24 : Igor Pavlov : Public domain */

#ifdef SHOW_STAT
static int ttt = 0;
#endif

static int verbosity = 0;

enum { 
        Fh_size = 6,    // file header size
        Ft_size = 20,   // file trailer size
                                        /*  0-3  CRC32 of the uncompressed data */
                                        /*  4-11 size of the uncompressed data */
                                        /* 12-19 member size including header and trailer */
};

typedef uint8_t File_trailer[Ft_size];

static inline void Ft_set_data_crc( File_trailer data, unsigned crc ) {
        for( int i = 0; i <= 3; ++i ) { data[i] = (uint8_t)crc; crc >>= 8; }
}

static inline void Ft_set_data_size( File_trailer data, unsigned long long sz ) {
        for( int i = 4; i <= 11; ++i ) { data[i] = (uint8_t)sz; sz >>= 8; }
}

static inline void Ft_set_member_size( File_trailer data, unsigned long long sz ) {
        for( int i = 12; i <= 19; ++i ) { data[i] = (uint8_t)sz; sz >>= 8; }
}

#define kNumTopBits 24
#define kTopValue ((uint32_t)1 << kNumTopBits)

#define kNumBitModelTotalBits 11
#define kBitModelTotal (1 << kNumBitModelTotalBits)
#define kNumMoveBits 5
#define kProbInitValue (kBitModelTotal >> 1)

#define kNumMoveReducingBits 4
#define kNumBitPriceShiftBits 4

#define kNumLogBits (9 + (int)sizeof(uint32_t) / 2)
#define kDicLogSizeMaxCompress ((kNumLogBits - 1) * 2 + 7)


static void LzmaEnc_FastPosInit(uint8_t *g_FastPos)
{
    int c = 2, slotFast;
    g_FastPos[0] = 0;
    g_FastPos[1] = 1;

    for (slotFast = 2; slotFast < kNumLogBits * 2; slotFast++)
    {
        uint32_t k = (1 << ((slotFast >> 1) - 1));
        uint32_t j;
        for (j = 0; j < k; j++, c++)
            g_FastPos[c] = (uint8_t)slotFast;
    }
}

#define BSR2_RET(pos, res) { uint32_t i = 6 + ((kNumLogBits - 1) & \
    (0 - (((((uint32_t)1 << (kNumLogBits + 6)) - 1) - pos) >> 31))); \
    res = p->g_FastPos[pos >> i] + (i * 2); }
/*
#define BSR2_RET(pos, res) { res = (pos < (1 << (kNumLogBits + 6))) ? \
    p->g_FastPos[pos >> 6] + 12 : \
    p->g_FastPos[pos >> (6 + kNumLogBits - 1)] + (6 + (kNumLogBits - 1)) * 2; }
*/

#define GetPosSlot1(pos) p->g_FastPos[pos]
#define GetPosSlot2(pos, res) { BSR2_RET(pos, res); }
#define GetPosSlot(pos, res) { if (pos < kNumFullDistances) res = p->g_FastPos[pos]; else BSR2_RET(pos, res); }


#define LZMA_NUM_REPS 4

typedef struct
{
    uint32_t price;

    State state;

    uint32_t posPrev2;
    uint32_t backPrev2;

    uint32_t posPrev;
    uint32_t backPrev;
    uint32_t backs[LZMA_NUM_REPS];

    bool prev1IsChar;
    bool prev2;
} COptimal;

#define kNumOpts (1 << 12)

#define kNumLenToPosStates 4
#define kNumPosSlotBits 6
#define kDicLogSizeMin 0
#define kDicLogSizeMax 32
#define kDistTableSizeMax (kDicLogSizeMax * 2)


#define kNumAlignBits 4
#define kAlignTableSize (1 << kNumAlignBits)
#define kAlignMask (kAlignTableSize - 1)

#define kStartPosModelIndex 4
#define kEndPosModelIndex 14
#define kNumPosModels (kEndPosModelIndex - kStartPosModelIndex)

#define kNumFullDistances (1 << (kEndPosModelIndex >> 1))

#define LZMA_LC_MAX 8
#define LZMA_LP_MAX 4
#define LZMA_PB_MAX 4

#define LZMA_NUM_PB_STATES_MAX (1 << LZMA_PB_MAX)


#define kLenNumLowBits 3
#define kLenNumLowSymbols (1 << kLenNumLowBits)
#define kLenNumMidBits 3
#define kLenNumMidSymbols (1 << kLenNumMidBits)
#define kLenNumHighBits 8
#define kLenNumHighSymbols (1 << kLenNumHighBits)

#define kLenNumSymbolsTotal (kLenNumLowSymbols + kLenNumMidSymbols + kLenNumHighSymbols)

#define LZMA_MATCH_LEN_MIN 2
#define LZMA_MATCH_LEN_MAX (LZMA_MATCH_LEN_MIN + kLenNumSymbolsTotal - 1)

#define kNumStates 12

typedef struct
{
    int choice;
    int choice2;
    int low[LZMA_NUM_PB_STATES_MAX << kLenNumLowBits];
    int mid[LZMA_NUM_PB_STATES_MAX << kLenNumMidBits];
    int high[kLenNumHighSymbols];
} CLenEnc;

typedef struct
{
    CLenEnc p;
    uint32_t prices[LZMA_NUM_PB_STATES_MAX][kLenNumSymbolsTotal];
    uint32_t tableSize;
    uint32_t counters[LZMA_NUM_PB_STATES_MAX];
} CLenPriceEnc;

typedef struct
{
    uint64_t low;
    uint64_t processed;
    uint8_t *bufBase;
    uint8_t *buf;
    uint8_t *bufLim;
    uint32_t range;
    uint32_t cacheSize;
    int outfd;
    int res;
    uint8_t cache;
} CRangeEnc;


typedef struct
{
    uint64_t nowPos64;
    int *litProbs;
    IMatchFinder matchFinder;
    CMatchFinder matchFinderBase;

    uint32_t optimumEndIndex;
    uint32_t optimumCurrentIndex;

    uint32_t longestMatchLength;
    uint32_t numPairs;
    uint32_t numAvail;
    COptimal opt[kNumOpts];

    uint8_t g_FastPos[1 << kNumLogBits];

    uint32_t ProbPrices[kBitModelTotal >> kNumMoveReducingBits];
    uint32_t matches[LZMA_MATCH_LEN_MAX * 2 + 2 + 1];
    uint32_t numFastBytes;
    uint32_t additionalOffset;
    uint32_t reps[LZMA_NUM_REPS];
    State state;

    uint32_t posSlotPrices[kNumLenToPosStates][kDistTableSizeMax];
    uint32_t distancesPrices[kNumLenToPosStates][kNumFullDistances];
    uint32_t alignPrices[kAlignTableSize];
    uint32_t alignPriceCount;

    uint32_t distTableSize;

    unsigned lc, lp, pb;
    unsigned lpMask, pbMask;

    int isMatch[kNumStates][LZMA_NUM_PB_STATES_MAX];
    int isRep[kNumStates];
    int isRepG0[kNumStates];
    int isRepG1[kNumStates];
    int isRepG2[kNumStates];
    int isRep0Long[kNumStates][LZMA_NUM_PB_STATES_MAX];

    int posSlotEncoder[kNumLenToPosStates][1 << kNumPosSlotBits];
    int posEncoders[kNumFullDistances - kEndPosModelIndex];
    int posAlignEncoder[1 << kNumAlignBits];

    CLenPriceEnc lenEnc;
    CLenPriceEnc repLenEnc;

    CRangeEnc rc;

    uint32_t matchPriceCount;

    int result;
    uint32_t dictSize;
    bool fastMode;
    bool finished;
} CLzmaEnc;


static const int kLiteralNextStates[kNumStates] = {0, 0, 0, 0, 1, 2, 3,  4,  5,  6,  4,  5};
static const int kMatchNextStates[kNumStates]   = {7, 7, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10};
static const int kRepNextStates[kNumStates]     = {8, 8, 8, 8, 8, 8, 8, 11, 11, 11, 11, 11};
static const int kShortRepNextStates[kNumStates]= {9, 9, 9, 9, 9, 9, 9, 11, 11, 11, 11, 11};

#define IsCharState(s) ((s) < 7)

#define GetLenToPosState(len) (((len) < kNumLenToPosStates + 1) ? (len) - 2 : kNumLenToPosStates - 1)

#define kInfinityPrice (1 << 30)

#define RC_BUF_SIZE (1 << 16)


static int RangeEnc_Init( CRangeEnc *p, const int outfd )
    {
    p->low = 0;
    p->processed = 0;
    p->range = 0xFFFFFFFF;
    p->cacheSize = 1;
    p->outfd = outfd;
    p->res = SZ_OK;
    p->cache = 0;
    p->buf = p->bufBase = (uint8_t *)LZMA_MALLOC( RC_BUF_SIZE );
    if( !p->bufBase ) return 0;
    p->bufLim = p->bufBase + RC_BUF_SIZE;
    return 1;
    }


static void RangeEnc_Free(CRangeEnc *p)
{
    LZMA_FREE(p->bufBase);
    p->bufBase = 0;
}


static void RangeEnc_FlushStream(CRangeEnc *p)
{
    int num;
    if (p->res != SZ_OK)
        return;
    num = p->buf - p->bufBase;
    if (num != writeblock(p->outfd, p->bufBase, num))
        p->res = SZ_ERROR_WRITE;
    p->processed += num;
    p->buf = p->bufBase;
}

static void RangeEnc_ShiftLow(CRangeEnc *p)
{
    if ((uint32_t)p->low < (uint32_t)0xFF000000 || (int)(p->low >> 32) != 0)
    {
        uint8_t temp = p->cache;
        do
        {
            uint8_t *buf = p->buf;
            *buf++ = (uint8_t)(temp + (uint8_t)(p->low >> 32));
            p->buf = buf;
            if (buf == p->bufLim)
                RangeEnc_FlushStream(p);
            temp = 0xFF;
        }
        while (--p->cacheSize != 0);
        p->cache = (uint8_t)((uint32_t)p->low >> 24);
    }
    p->cacheSize++;
    p->low = (uint32_t)p->low << 8;
}

static void RangeEnc_FlushData(CRangeEnc *p)
{
    int i;
    for (i = 0; i < 5; i++)
        RangeEnc_ShiftLow(p);
}

static void RangeEnc_EncodeDirectBits(CRangeEnc *p, uint32_t value, int numBits)
{
    do
    {
        p->range >>= 1;
        p->low += p->range & (0 - ((value >> --numBits) & 1));
        if (p->range < kTopValue)
        {
            p->range <<= 8;
            RangeEnc_ShiftLow(p);
        }
    }
    while (numBits != 0);
}

static void RangeEnc_EncodeBit(CRangeEnc *p, int *prob, uint32_t symbol)
{
    uint32_t ttt = *prob;
    uint32_t newBound = (p->range >> kNumBitModelTotalBits) * ttt;
    if (symbol == 0)
    {
        p->range = newBound;
        ttt += (kBitModelTotal - ttt) >> kNumMoveBits;
    }
    else
    {
        p->low += newBound;
        p->range -= newBound;
        ttt -= ttt >> kNumMoveBits;
    }
    *prob = (int)ttt;
    if (p->range < kTopValue)
    {
        p->range <<= 8;
        RangeEnc_ShiftLow(p);
    }
}

static void LitEnc_Encode(CRangeEnc *p, int *probs, uint32_t symbol)
{
    symbol |= 0x100;
    do
    {
        RangeEnc_EncodeBit(p, probs + (symbol >> 8), (symbol >> 7) & 1);
        symbol <<= 1;
    }
    while (symbol < 0x10000);
}

static void LitEnc_EncodeMatched(CRangeEnc *p, int *probs, uint32_t symbol, uint32_t matchByte)
{
    uint32_t offs = 0x100;
    symbol |= 0x100;
    do
    {
        matchByte <<= 1;
        RangeEnc_EncodeBit(p, probs + (offs + (matchByte & offs) + (symbol >> 8)), (symbol >> 7) & 1);
        symbol <<= 1;
        offs &= ~(matchByte ^ symbol);
    }
    while (symbol < 0x10000);
}

static void LzmaEnc_InitPriceTables(uint32_t *ProbPrices)
{
    uint32_t i;
    for (i = (1 << kNumMoveReducingBits) / 2; i < kBitModelTotal; i += (1 << kNumMoveReducingBits))
    {
        const int kCyclesBits = kNumBitPriceShiftBits;
        uint32_t w = i;
        uint32_t bitCount = 0;
        int j;
        for (j = 0; j < kCyclesBits; j++)
        {
            w = w * w;
            bitCount <<= 1;
            while (w >= ((uint32_t)1 << 16))
            {
                w >>= 1;
                bitCount++;
            }
        }
        ProbPrices[i >> kNumMoveReducingBits] = ((kNumBitModelTotalBits << kCyclesBits) - 15 - bitCount);
    }
}


#define GET_PRICE(prob, symbol) \
    p->ProbPrices[((prob) ^ (((-(int)(symbol))) & (kBitModelTotal - 1))) >> kNumMoveReducingBits];

#define GET_PRICEa(prob, symbol) \
    ProbPrices[((prob) ^ ((-((int)(symbol))) & (kBitModelTotal - 1))) >> kNumMoveReducingBits];

#define GET_PRICE_0(prob) p->ProbPrices[(prob) >> kNumMoveReducingBits]
#define GET_PRICE_1(prob) p->ProbPrices[((prob) ^ (kBitModelTotal - 1)) >> kNumMoveReducingBits]

#define GET_PRICE_0a(prob) ProbPrices[(prob) >> kNumMoveReducingBits]
#define GET_PRICE_1a(prob) ProbPrices[((prob) ^ (kBitModelTotal - 1)) >> kNumMoveReducingBits]

static uint32_t LitEnc_GetPrice(const int *probs, uint32_t symbol, uint32_t *ProbPrices)
{
    uint32_t price = 0;
    symbol |= 0x100;
    do
    {
        price += GET_PRICEa(probs[symbol >> 8], (symbol >> 7) & 1);
        symbol <<= 1;
    }
    while (symbol < 0x10000);
    return price;
}

static uint32_t LitEnc_GetPriceMatched(const int *probs, uint32_t symbol, uint32_t matchByte, uint32_t *ProbPrices)
{
    uint32_t price = 0;
    uint32_t offs = 0x100;
    symbol |= 0x100;
    do
    {
        matchByte <<= 1;
        price += GET_PRICEa(probs[offs + (matchByte & offs) + (symbol >> 8)], (symbol >> 7) & 1);
        symbol <<= 1;
        offs &= ~(matchByte ^ symbol);
    }
    while (symbol < 0x10000);
    return price;
}


static void RcTree_Encode(CRangeEnc *rc, int *probs, int numBitLevels, uint32_t symbol)
{
    uint32_t m = 1;
    int i;
    for (i = numBitLevels; i != 0;)
    {
        uint32_t bit;
        i--;
        bit = (symbol >> i) & 1;
        RangeEnc_EncodeBit(rc, probs + m, bit);
        m = (m << 1) | bit;
    }
}

static void RcTree_ReverseEncode(CRangeEnc *rc, int *probs, int numBitLevels, uint32_t symbol)
{
    uint32_t m = 1;
    int i;
    for (i = 0; i < numBitLevels; i++)
    {
        uint32_t bit = symbol & 1;
        RangeEnc_EncodeBit(rc, probs + m, bit);
        m = (m << 1) | bit;
        symbol >>= 1;
    }
}

static uint32_t RcTree_GetPrice(const int *probs, int numBitLevels, uint32_t symbol, uint32_t *ProbPrices)
{
    uint32_t price = 0;
    symbol |= (1 << numBitLevels);
    while (symbol != 1)
    {
        price += GET_PRICEa(probs[symbol >> 1], symbol & 1);
        symbol >>= 1;
    }
    return price;
}

static uint32_t RcTree_ReverseGetPrice(const int *probs, int numBitLevels, uint32_t symbol, uint32_t *ProbPrices)
{
    uint32_t price = 0;
    uint32_t m = 1;
    int i;
    for (i = numBitLevels; i != 0; i--)
    {
        uint32_t bit = symbol & 1;
        symbol >>= 1;
        price += GET_PRICEa(probs[m], bit);
        m = (m << 1) | bit;
    }
    return price;
}


static void LenEnc_Init(CLenEnc *p)
{
    unsigned i;
    p->choice = p->choice2 = kProbInitValue;
    for (i = 0; i < (LZMA_NUM_PB_STATES_MAX << kLenNumLowBits); i++)
        p->low[i] = kProbInitValue;
    for (i = 0; i < (LZMA_NUM_PB_STATES_MAX << kLenNumMidBits); i++)
        p->mid[i] = kProbInitValue;
    for (i = 0; i < kLenNumHighSymbols; i++)
        p->high[i] = kProbInitValue;
}

static void LenEnc_Encode(CLenEnc *p, CRangeEnc *rc, uint32_t symbol, uint32_t posState)
{
    if (symbol < kLenNumLowSymbols)
    {
        RangeEnc_EncodeBit(rc, &p->choice, 0);
        RcTree_Encode(rc, p->low + (posState << kLenNumLowBits), kLenNumLowBits, symbol);
    }
    else
    {
        RangeEnc_EncodeBit(rc, &p->choice, 1);
        if (symbol < kLenNumLowSymbols + kLenNumMidSymbols)
        {
            RangeEnc_EncodeBit(rc, &p->choice2, 0);
            RcTree_Encode(rc, p->mid + (posState << kLenNumMidBits), kLenNumMidBits, symbol - kLenNumLowSymbols);
        }
        else
        {
            RangeEnc_EncodeBit(rc, &p->choice2, 1);
            RcTree_Encode(rc, p->high, kLenNumHighBits, symbol - kLenNumLowSymbols - kLenNumMidSymbols);
        }
    }
}

static void LenEnc_SetPrices(CLenEnc *p, uint32_t posState, uint32_t numSymbols, uint32_t *prices, uint32_t *ProbPrices)
{
    uint32_t a0 = GET_PRICE_0a(p->choice);
    uint32_t a1 = GET_PRICE_1a(p->choice);
    uint32_t b0 = a1 + GET_PRICE_0a(p->choice2);
    uint32_t b1 = a1 + GET_PRICE_1a(p->choice2);
    uint32_t i = 0;
    for (i = 0; i < kLenNumLowSymbols; i++)
    {
        if (i >= numSymbols)
            return;
        prices[i] = a0 + RcTree_GetPrice(p->low + (posState << kLenNumLowBits), kLenNumLowBits, i, ProbPrices);
    }
    for (; i < kLenNumLowSymbols + kLenNumMidSymbols; i++)
    {
        if (i >= numSymbols)
            return;
        prices[i] = b0 + RcTree_GetPrice(p->mid + (posState << kLenNumMidBits), kLenNumMidBits, i - kLenNumLowSymbols, ProbPrices);
    }
    for (; i < numSymbols; i++)
        prices[i] = b1 + RcTree_GetPrice(p->high, kLenNumHighBits, i - kLenNumLowSymbols - kLenNumMidSymbols, ProbPrices);
}

static void LenPriceEnc_UpdateTable(CLenPriceEnc *p, uint32_t posState, uint32_t *ProbPrices)
{
    LenEnc_SetPrices(&p->p, posState, p->tableSize, p->prices[posState], ProbPrices);
    p->counters[posState] = p->tableSize;
}

static void LenPriceEnc_UpdateTables(CLenPriceEnc *p, uint32_t numPosStates, uint32_t *ProbPrices)
{
    uint32_t posState;
    for (posState = 0; posState < numPosStates; posState++)
        LenPriceEnc_UpdateTable(p, posState, ProbPrices);
}

static void LenEnc_Encode2(CLenPriceEnc *p, CRangeEnc *rc, uint32_t symbol, uint32_t posState, bool updatePrice, uint32_t *ProbPrices)
{
    LenEnc_Encode(&p->p, rc, symbol, posState);
    if (updatePrice)
        if (--p->counters[posState] == 0)
            LenPriceEnc_UpdateTable(p, posState, ProbPrices);
}


static void MovePos(CLzmaEnc *p, uint32_t num)
{
    #ifdef SHOW_STAT
    ttt += num;
    printf("\n MovePos %d", num);
    #endif
    if (num != 0)
    {
        p->additionalOffset += num;
        p->matchFinder.Skip(&p->matchFinderBase, num);
    }
}

static uint32_t ReadMatchDistances(CLzmaEnc *p, uint32_t *numDistancePairsRes)
{
    uint32_t lenRes = 0, numPairs;
    p->numAvail = Mf_GetNumAvailableBytes(&p->matchFinderBase);
    numPairs = p->matchFinder.GetMatches(&p->matchFinderBase, p->matches);
    #ifdef SHOW_STAT
    printf("\n i = %d numPairs = %d    ", ttt, numPairs / 2);
    ttt++;
    {
        uint32_t i;
        for (i = 0; i < numPairs; i += 2)
            printf("%2d %6d   | ", p->matches[i], p->matches[i + 1]);
    }
    #endif
    if (numPairs > 0)
    {
        lenRes = p->matches[numPairs - 2];
        if (lenRes == p->numFastBytes)
        {
            const uint8_t *pby = Mf_GetPointerToCurrentPos(&p->matchFinderBase) - 1;
            uint32_t distance = p->matches[numPairs - 1] + 1;
            uint32_t numAvail = p->numAvail;
            if (numAvail > LZMA_MATCH_LEN_MAX)
                numAvail = LZMA_MATCH_LEN_MAX;
            {
                const uint8_t *pby2 = pby - distance;
                for (; lenRes < numAvail && pby[lenRes] == pby2[lenRes]; lenRes++) ;
            }
        }
    }
    p->additionalOffset++;
    *numDistancePairsRes = numPairs;
    return lenRes;
}


#define MakeAsChar(p) (p)->backPrev = (uint32_t)(-1); (p)->prev1IsChar = false;
#define MakeAsShortRep(p) (p)->backPrev = 0; (p)->prev1IsChar = false;
#define IsShortRep(p) ((p)->backPrev == 0)

static uint32_t GetRepLen1Price(CLzmaEnc *p, State state, uint32_t posState)
{
    return
        GET_PRICE_0(p->isRepG0[state]) +
        GET_PRICE_0(p->isRep0Long[state][posState]);
}

static uint32_t GetPureRepPrice(CLzmaEnc *p, uint32_t repIndex, State state, uint32_t posState)
{
    uint32_t price;
    if (repIndex == 0)
    {
        price = GET_PRICE_0(p->isRepG0[state]);
        price += GET_PRICE_1(p->isRep0Long[state][posState]);
    }
    else
    {
        price = GET_PRICE_1(p->isRepG0[state]);
        if (repIndex == 1)
            price += GET_PRICE_0(p->isRepG1[state]);
        else
        {
            price += GET_PRICE_1(p->isRepG1[state]);
            price += GET_PRICE(p->isRepG2[state], repIndex - 2);
        }
    }
    return price;
}

static uint32_t GetRepPrice(CLzmaEnc *p, uint32_t repIndex, uint32_t len, State state, uint32_t posState)
{
    return p->repLenEnc.prices[posState][len - LZMA_MATCH_LEN_MIN] +
        GetPureRepPrice(p, repIndex, state, posState);
}

static uint32_t Backward(CLzmaEnc *p, uint32_t *backRes, uint32_t cur)
{
    uint32_t posMem = p->opt[cur].posPrev;
    uint32_t backMem = p->opt[cur].backPrev;
    p->optimumEndIndex = cur;
    do
    {
        if (p->opt[cur].prev1IsChar)
        {
            MakeAsChar(&p->opt[posMem])
            p->opt[posMem].posPrev = posMem - 1;
            if (p->opt[cur].prev2)
            {
                p->opt[posMem - 1].prev1IsChar = false;
                p->opt[posMem - 1].posPrev = p->opt[cur].posPrev2;
                p->opt[posMem - 1].backPrev = p->opt[cur].backPrev2;
            }
        }
        {
            uint32_t posPrev = posMem;
            uint32_t backCur = backMem;

            backMem = p->opt[posPrev].backPrev;
            posMem = p->opt[posPrev].posPrev;

            p->opt[posPrev].backPrev = backCur;
            p->opt[posPrev].posPrev = cur;
            cur = posPrev;
        }
    }
    while (cur != 0);
    *backRes = p->opt[0].backPrev;
    p->optimumCurrentIndex  = p->opt[0].posPrev;
    return p->optimumCurrentIndex;
}

#define LIT_PROBS(pos, prevByte) (p->litProbs + ((((pos) & p->lpMask) << p->lc) + ((prevByte) >> (8 - p->lc))) * 0x300)

static uint32_t GetOptimum(CLzmaEnc *p, uint32_t position, uint32_t *backRes)
{
    uint32_t numAvail, mainLen, numPairs, repMaxIndex, i, posState, lenEnd, len, cur;
    uint32_t matchPrice, repMatchPrice, normalMatchPrice;
    uint32_t reps[LZMA_NUM_REPS], repLens[LZMA_NUM_REPS];
    uint32_t *matches;
    const uint8_t *data;
    uint8_t curByte, matchByte;
    if (p->optimumEndIndex != p->optimumCurrentIndex)
    {
        const COptimal *opt = &p->opt[p->optimumCurrentIndex];
        uint32_t lenRes = opt->posPrev - p->optimumCurrentIndex;
        *backRes = opt->backPrev;
        p->optimumCurrentIndex = opt->posPrev;
        return lenRes;
    }
    p->optimumCurrentIndex = p->optimumEndIndex = 0;

    if (p->additionalOffset == 0)
        mainLen = ReadMatchDistances(p, &numPairs);
    else
    {
        mainLen = p->longestMatchLength;
        numPairs = p->numPairs;
    }

    numAvail = p->numAvail;
    if (numAvail < 2)
    {
        *backRes = (uint32_t)(-1);
        return 1;
    }
    if (numAvail > LZMA_MATCH_LEN_MAX)
        numAvail = LZMA_MATCH_LEN_MAX;

    data = Mf_GetPointerToCurrentPos(&p->matchFinderBase) - 1;
    repMaxIndex = 0;
    for (i = 0; i < LZMA_NUM_REPS; i++)
    {
        uint32_t lenTest;
        const uint8_t *data2;
        reps[i] = p->reps[i];
        data2 = data - (reps[i] + 1);
        if (data[0] != data2[0] || data[1] != data2[1])
        {
            repLens[i] = 0;
            continue;
        }
        for (lenTest = 2; lenTest < numAvail && data[lenTest] == data2[lenTest]; lenTest++) ;
        repLens[i] = lenTest;
        if (lenTest > repLens[repMaxIndex])
            repMaxIndex = i;
    }
    if (repLens[repMaxIndex] >= p->numFastBytes)
    {
        uint32_t lenRes;
        *backRes = repMaxIndex;
        lenRes = repLens[repMaxIndex];
        MovePos(p, lenRes - 1);
        return lenRes;
    }

    matches = p->matches;
    if (mainLen >= p->numFastBytes)
    {
        *backRes = matches[numPairs - 1] + LZMA_NUM_REPS;
        MovePos(p, mainLen - 1);
        return mainLen;
    }
    curByte = *data;
    matchByte = *(data - (reps[0] + 1));

    if (mainLen < 2 && curByte != matchByte && repLens[repMaxIndex] < 2)
    {
        *backRes = (uint32_t)-1;
        return 1;
    }

    p->opt[0].state = p->state;

    posState = (position & p->pbMask);

    {
        const int *probs = LIT_PROBS(position, *(data - 1));
        p->opt[1].price = GET_PRICE_0(p->isMatch[p->state][posState]) +
                (!IsCharState(p->state) ?
                    LitEnc_GetPriceMatched(probs, curByte, matchByte, p->ProbPrices) :
                    LitEnc_GetPrice(probs, curByte, p->ProbPrices));
    }

    MakeAsChar(&p->opt[1]);

    matchPrice = GET_PRICE_1(p->isMatch[p->state][posState]);
    repMatchPrice = matchPrice + GET_PRICE_1(p->isRep[p->state]);

    if (matchByte == curByte)
    {
        uint32_t shortRepPrice = repMatchPrice + GetRepLen1Price(p, p->state, posState);
        if (shortRepPrice < p->opt[1].price)
        {
            p->opt[1].price = shortRepPrice;
            MakeAsShortRep(&p->opt[1]);
        }
    }
    lenEnd = ((mainLen >= repLens[repMaxIndex]) ? mainLen : repLens[repMaxIndex]);

    if (lenEnd < 2)
    {
        *backRes = p->opt[1].backPrev;
        return 1;
    }

    p->opt[1].posPrev = 0;
    for (i = 0; i < LZMA_NUM_REPS; i++)
        p->opt[0].backs[i] = reps[i];

    len = lenEnd;
    do
        p->opt[len--].price = kInfinityPrice;
    while (len >= 2);

    for (i = 0; i < LZMA_NUM_REPS; i++)
    {
        uint32_t repLen = repLens[i];
        uint32_t price;
        if (repLen < 2)
            continue;
        price = repMatchPrice + GetPureRepPrice(p, i, p->state, posState);
        do
        {
            uint32_t curAndLenPrice = price + p->repLenEnc.prices[posState][repLen - 2];
            COptimal *opt = &p->opt[repLen];
            if (curAndLenPrice < opt->price)
            {
                opt->price = curAndLenPrice;
                opt->posPrev = 0;
                opt->backPrev = i;
                opt->prev1IsChar = false;
            }
        }
        while (--repLen >= 2);
    }

    normalMatchPrice = matchPrice + GET_PRICE_0(p->isRep[p->state]);

    len = ((repLens[0] >= 2) ? repLens[0] + 1 : 2);
    if (len <= mainLen)
    {
        uint32_t offs = 0;
        while (len > matches[offs])
            offs += 2;
        for (; ; len++)
        {
            COptimal *opt;
            uint32_t distance = matches[offs + 1];

            uint32_t curAndLenPrice = normalMatchPrice + p->lenEnc.prices[posState][len - LZMA_MATCH_LEN_MIN];
            uint32_t lenToPosState = GetLenToPosState(len);
            if (distance < kNumFullDistances)
                curAndLenPrice += p->distancesPrices[lenToPosState][distance];
            else
            {
                uint32_t slot;
                GetPosSlot2(distance, slot);
                curAndLenPrice += p->alignPrices[distance & kAlignMask] + p->posSlotPrices[lenToPosState][slot];
            }
            opt = &p->opt[len];
            if (curAndLenPrice < opt->price)
            {
                opt->price = curAndLenPrice;
                opt->posPrev = 0;
                opt->backPrev = distance + LZMA_NUM_REPS;
                opt->prev1IsChar = false;
            }
            if (len == matches[offs])
            {
                offs += 2;
                if (offs == numPairs)
                    break;
            }
        }
    }

    cur = 0;

        #ifdef SHOW_STAT2
        if (position >= 0)
        {
            unsigned i;
            printf("\n pos = %4X", position);
            for (i = cur; i <= lenEnd; i++)
            printf("\nprice[%4X] = %d", position - cur + i, p->opt[i].price);
        }
        #endif

    for (;;)
    {
        uint32_t numAvailFull, newLen, numPairs, posPrev, state, posState, startLen;
        uint32_t curPrice, curAnd1Price, matchPrice, repMatchPrice;
        bool nextIsChar;
        uint8_t curByte, matchByte;
        const uint8_t *data;
        COptimal *curOpt;
        COptimal *nextOpt;

        cur++;
        if (cur == lenEnd)
            return Backward(p, backRes, cur);

        newLen = ReadMatchDistances(p, &numPairs);
        if (newLen >= p->numFastBytes)
        {
            p->numPairs = numPairs;
            p->longestMatchLength = newLen;
            return Backward(p, backRes, cur);
        }
        position++;
        curOpt = &p->opt[cur];
        posPrev = curOpt->posPrev;
        if (curOpt->prev1IsChar)
        {
            posPrev--;
            if (curOpt->prev2)
            {
                state = p->opt[curOpt->posPrev2].state;
                if (curOpt->backPrev2 < LZMA_NUM_REPS)
                    state = kRepNextStates[state];
                else
                    state = kMatchNextStates[state];
            }
            else
                state = p->opt[posPrev].state;
            state = kLiteralNextStates[state];
        }
        else
            state = p->opt[posPrev].state;
        if (posPrev == cur - 1)
        {
            if (IsShortRep(curOpt))
                state = kShortRepNextStates[state];
            else
                state = kLiteralNextStates[state];
        }
        else
        {
            uint32_t pos;
            const COptimal *prevOpt;
            if (curOpt->prev1IsChar && curOpt->prev2)
            {
                posPrev = curOpt->posPrev2;
                pos = curOpt->backPrev2;
                state = kRepNextStates[state];
            }
            else
            {
                pos = curOpt->backPrev;
                if (pos < LZMA_NUM_REPS)
                    state = kRepNextStates[state];
                else
                    state = kMatchNextStates[state];
            }
            prevOpt = &p->opt[posPrev];
            if (pos < LZMA_NUM_REPS)
            {
                uint32_t i;
                reps[0] = prevOpt->backs[pos];
                for (i = 1; i <= pos; i++)
                    reps[i] = prevOpt->backs[i - 1];
                for (; i < LZMA_NUM_REPS; i++)
                    reps[i] = prevOpt->backs[i];
            }
            else
            {
                uint32_t i;
                reps[0] = (pos - LZMA_NUM_REPS);
                for (i = 1; i < LZMA_NUM_REPS; i++)
                    reps[i] = prevOpt->backs[i - 1];
            }
        }
        curOpt->state = state;

        curOpt->backs[0] = reps[0];
        curOpt->backs[1] = reps[1];
        curOpt->backs[2] = reps[2];
        curOpt->backs[3] = reps[3];

        curPrice = curOpt->price;
        nextIsChar = false;
        data = Mf_GetPointerToCurrentPos(&p->matchFinderBase) - 1;
        curByte = *data;
        matchByte = *(data - (reps[0] + 1));

        posState = (position & p->pbMask);

        curAnd1Price = curPrice + GET_PRICE_0(p->isMatch[state][posState]);
        {
            const int *probs = LIT_PROBS(position, *(data - 1));
            curAnd1Price +=
                (!IsCharState(state) ?
                    LitEnc_GetPriceMatched(probs, curByte, matchByte, p->ProbPrices) :
                    LitEnc_GetPrice(probs, curByte, p->ProbPrices));
        }

        nextOpt = &p->opt[cur + 1];

        if (curAnd1Price < nextOpt->price)
        {
            nextOpt->price = curAnd1Price;
            nextOpt->posPrev = cur;
            MakeAsChar(nextOpt);
            nextIsChar = true;
        }

        matchPrice = curPrice + GET_PRICE_1(p->isMatch[state][posState]);
        repMatchPrice = matchPrice + GET_PRICE_1(p->isRep[state]);

        if (matchByte == curByte && !(nextOpt->posPrev < cur && nextOpt->backPrev == 0))
        {
            uint32_t shortRepPrice = repMatchPrice + GetRepLen1Price(p, state, posState);
            if (shortRepPrice <= nextOpt->price)
            {
                nextOpt->price = shortRepPrice;
                nextOpt->posPrev = cur;
                MakeAsShortRep(nextOpt);
                nextIsChar = true;
            }
        }
        numAvailFull = p->numAvail;
        {
            uint32_t temp = kNumOpts - 1 - cur;
            if (temp < numAvailFull)
                numAvailFull = temp;
        }

        if (numAvailFull < 2)
            continue;
        numAvail = (numAvailFull <= p->numFastBytes ? numAvailFull : p->numFastBytes);

        if (!nextIsChar && matchByte != curByte) /* speed optimization */
        {
            /* try Literal + rep0 */
            uint32_t temp;
            uint32_t lenTest2;
            const uint8_t *data2 = data - (reps[0] + 1);
            uint32_t limit = p->numFastBytes + 1;
            if (limit > numAvailFull)
                limit = numAvailFull;

            for (temp = 1; temp < limit && data[temp] == data2[temp]; temp++) ;
            lenTest2 = temp - 1;
            if (lenTest2 >= 2)
            {
                State state2 = kLiteralNextStates[state];
                uint32_t posStateNext = (position + 1) & p->pbMask;
                uint32_t nextRepMatchPrice = curAnd1Price +
                        GET_PRICE_1(p->isMatch[state2][posStateNext]) +
                        GET_PRICE_1(p->isRep[state2]);
                /* for (; lenTest2 >= 2; lenTest2--) */
                {
                    uint32_t curAndLenPrice;
                    COptimal *opt;
                    uint32_t offset = cur + 1 + lenTest2;
                    while (lenEnd < offset)
                        p->opt[++lenEnd].price = kInfinityPrice;
                    curAndLenPrice = nextRepMatchPrice + GetRepPrice(p, 0, lenTest2, state2, posStateNext);
                    opt = &p->opt[offset];
                    if (curAndLenPrice < opt->price)
                    {
                        opt->price = curAndLenPrice;
                        opt->posPrev = cur + 1;
                        opt->backPrev = 0;
                        opt->prev1IsChar = true;
                        opt->prev2 = false;
                    }
                }
            }
        }

        startLen = 2; /* speed optimization */
        {
        uint32_t repIndex;
        for (repIndex = 0; repIndex < LZMA_NUM_REPS; repIndex++)
        {
            uint32_t lenTest;
            uint32_t lenTestTemp;
            uint32_t price;
            const uint8_t *data2 = data - (reps[repIndex] + 1);
            if (data[0] != data2[0] || data[1] != data2[1])
                continue;
            for (lenTest = 2; lenTest < numAvail && data[lenTest] == data2[lenTest]; lenTest++) ;
            while (lenEnd < cur + lenTest)
                p->opt[++lenEnd].price = kInfinityPrice;
            lenTestTemp = lenTest;
            price = repMatchPrice + GetPureRepPrice(p, repIndex, state, posState);
            do
            {
                uint32_t curAndLenPrice = price + p->repLenEnc.prices[posState][lenTest - 2];
                COptimal *opt = &p->opt[cur + lenTest];
                if (curAndLenPrice < opt->price)
                {
                    opt->price = curAndLenPrice;
                    opt->posPrev = cur;
                    opt->backPrev = repIndex;
                    opt->prev1IsChar = false;
                }
            }
            while (--lenTest >= 2);
            lenTest = lenTestTemp;

            if (repIndex == 0)
                startLen = lenTest + 1;

            /* if (_maxMode) */
                {
                    uint32_t lenTest2 = lenTest + 1;
                    uint32_t limit = lenTest2 + p->numFastBytes;
                    uint32_t nextRepMatchPrice;
                    if (limit > numAvailFull)
                        limit = numAvailFull;
                    for (; lenTest2 < limit && data[lenTest2] == data2[lenTest2]; lenTest2++) ;
                    lenTest2 -= lenTest + 1;
                    if (lenTest2 >= 2)
                    {
                        State state2 = kRepNextStates[state];
                        uint32_t posStateNext = (position + lenTest) & p->pbMask;
                        uint32_t curAndLenCharPrice =
                                price + p->repLenEnc.prices[posState][lenTest - 2] +
                                GET_PRICE_0(p->isMatch[state2][posStateNext]) +
                                LitEnc_GetPriceMatched(LIT_PROBS(position + lenTest, data[lenTest - 1]),
                                        data[lenTest], data2[lenTest], p->ProbPrices);
                        state2 = kLiteralNextStates[state2];
                        posStateNext = (position + lenTest + 1) & p->pbMask;
                        nextRepMatchPrice = curAndLenCharPrice +
                                GET_PRICE_1(p->isMatch[state2][posStateNext]) +
                                GET_PRICE_1(p->isRep[state2]);

                        /* for (; lenTest2 >= 2; lenTest2--) */
                        {
                            uint32_t curAndLenPrice;
                            COptimal *opt;
                            uint32_t offset = cur + lenTest + 1 + lenTest2;
                            while (lenEnd < offset)
                                p->opt[++lenEnd].price = kInfinityPrice;
                            curAndLenPrice = nextRepMatchPrice + GetRepPrice(p, 0, lenTest2, state2, posStateNext);
                            opt = &p->opt[offset];
                            if (curAndLenPrice < opt->price)
                            {
                                opt->price = curAndLenPrice;
                                opt->posPrev = cur + lenTest + 1;
                                opt->backPrev = 0;
                                opt->prev1IsChar = true;
                                opt->prev2 = true;
                                opt->posPrev2 = cur;
                                opt->backPrev2 = repIndex;
                            }
                        }
                    }
                }
        }
        }
        /* for (uint32_t lenTest = 2; lenTest <= newLen; lenTest++) */
        if (newLen > numAvail)
        {
            newLen = numAvail;
            for (numPairs = 0; newLen > matches[numPairs]; numPairs += 2) ;
            matches[numPairs] = newLen;
            numPairs += 2;
        }
        if (newLen >= startLen)
        {
            uint32_t normalMatchPrice = matchPrice + GET_PRICE_0(p->isRep[state]);
            uint32_t offs, curBack, posSlot;
            uint32_t lenTest;
            while (lenEnd < cur + newLen)
                p->opt[++lenEnd].price = kInfinityPrice;

            offs = 0;
            while (startLen > matches[offs])
                offs += 2;
            curBack = matches[offs + 1];
            GetPosSlot2(curBack, posSlot);
            for (lenTest = /*2*/ startLen; ; lenTest++)
            {
                uint32_t curAndLenPrice = normalMatchPrice + p->lenEnc.prices[posState][lenTest - LZMA_MATCH_LEN_MIN];
                uint32_t lenToPosState = GetLenToPosState(lenTest);
                COptimal *opt;
                if (curBack < kNumFullDistances)
                    curAndLenPrice += p->distancesPrices[lenToPosState][curBack];
                else
                    curAndLenPrice += p->posSlotPrices[lenToPosState][posSlot] + p->alignPrices[curBack & kAlignMask];

                opt = &p->opt[cur + lenTest];
                if (curAndLenPrice < opt->price)
                {
                    opt->price = curAndLenPrice;
                    opt->posPrev = cur;
                    opt->backPrev = curBack + LZMA_NUM_REPS;
                    opt->prev1IsChar = false;
                }

                if (/*_maxMode && */lenTest == matches[offs])
                {
                    /* Try Match + Literal + Rep0 */
                    const uint8_t *data2 = data - (curBack + 1);
                    uint32_t lenTest2 = lenTest + 1;
                    uint32_t limit = lenTest2 + p->numFastBytes;
                    uint32_t nextRepMatchPrice;
                    if (limit > numAvailFull)
                        limit = numAvailFull;
                    for (; lenTest2 < limit && data[lenTest2] == data2[lenTest2]; lenTest2++) ;
                    lenTest2 -= lenTest + 1;
                    if (lenTest2 >= 2)
                    {
                        State state2 = kMatchNextStates[state];
                        uint32_t posStateNext = (position + lenTest) & p->pbMask;
                        uint32_t curAndLenCharPrice = curAndLenPrice +
                                GET_PRICE_0(p->isMatch[state2][posStateNext]) +
                                LitEnc_GetPriceMatched(LIT_PROBS(position + lenTest, data[lenTest - 1]),
                                        data[lenTest], data2[lenTest], p->ProbPrices);
                        state2 = kLiteralNextStates[state2];
                        posStateNext = (posStateNext + 1) & p->pbMask;
                        nextRepMatchPrice = curAndLenCharPrice +
                                GET_PRICE_1(p->isMatch[state2][posStateNext]) +
                                GET_PRICE_1(p->isRep[state2]);

                        /* for (; lenTest2 >= 2; lenTest2--) */
                        {
                            uint32_t offset = cur + lenTest + 1 + lenTest2;
                            uint32_t curAndLenPrice;
                            COptimal *opt;
                            while (lenEnd < offset)
                                p->opt[++lenEnd].price = kInfinityPrice;
                            curAndLenPrice = nextRepMatchPrice + GetRepPrice(p, 0, lenTest2, state2, posStateNext);
                            opt = &p->opt[offset];
                            if (curAndLenPrice < opt->price)
                            {
                                opt->price = curAndLenPrice;
                                opt->posPrev = cur + lenTest + 1;
                                opt->backPrev = 0;
                                opt->prev1IsChar = true;
                                opt->prev2 = true;
                                opt->posPrev2 = cur;
                                opt->backPrev2 = curBack + LZMA_NUM_REPS;
                            }
                        }
                    }
                    offs += 2;
                    if (offs == numPairs)
                        break;
                    curBack = matches[offs + 1];
                    if (curBack >= kNumFullDistances)
                        GetPosSlot2(curBack, posSlot);
                }
            }
        }
    }
}

#define ChangePair(smallDist, bigDist) (((bigDist) >> 7) > (smallDist))

static uint32_t GetOptimumFast(CLzmaEnc *p, uint32_t *backRes)
{
    uint32_t numAvail, mainLen, mainDist, numPairs, repIndex, repLen, i;
    const uint8_t *data;
    const uint32_t *matches;

    if (p->additionalOffset == 0)
        mainLen = ReadMatchDistances(p, &numPairs);
    else
    {
        mainLen = p->longestMatchLength;
        numPairs = p->numPairs;
    }

    numAvail = p->numAvail;
    *backRes = (uint32_t)-1;
    if (numAvail < 2)
        return 1;
    if (numAvail > LZMA_MATCH_LEN_MAX)
        numAvail = LZMA_MATCH_LEN_MAX;
    data = Mf_GetPointerToCurrentPos(&p->matchFinderBase) - 1;

    repLen = repIndex = 0;
    for (i = 0; i < LZMA_NUM_REPS; i++)
    {
        uint32_t len;
        const uint8_t *data2 = data - (p->reps[i] + 1);
        if (data[0] != data2[0] || data[1] != data2[1])
            continue;
        for (len = 2; len < numAvail && data[len] == data2[len]; len++) ;
        if (len >= p->numFastBytes)
        {
            *backRes = i;
            MovePos(p, len - 1);
            return len;
        }
        if (len > repLen)
        {
            repIndex = i;
            repLen = len;
        }
    }

    matches = p->matches;
    if (mainLen >= p->numFastBytes)
    {
        *backRes = matches[numPairs - 1] + LZMA_NUM_REPS;
        MovePos(p, mainLen - 1);
        return mainLen;
    }

    mainDist = 0; /* for GCC */
    if (mainLen >= 2)
    {
        mainDist = matches[numPairs - 1];
        while (numPairs > 2 && mainLen == matches[numPairs - 4] + 1)
        {
            if (!ChangePair(matches[numPairs - 3], mainDist))
                break;
            numPairs -= 2;
            mainLen = matches[numPairs - 2];
            mainDist = matches[numPairs - 1];
        }
        if (mainLen == 2 && mainDist >= 0x80)
            mainLen = 1;
    }

    if (repLen >= 2 && (
                (repLen + 1 >= mainLen) ||
                (repLen + 2 >= mainLen && mainDist >= (1 << 9)) ||
                (repLen + 3 >= mainLen && mainDist >= (1 << 15))))
    {
        *backRes = repIndex;
        MovePos(p, repLen - 1);
        return repLen;
    }

    if (mainLen < 2 || numAvail <= 2)
        return 1;

    p->longestMatchLength = ReadMatchDistances(p, &p->numPairs);
    if (p->longestMatchLength >= 2)
    {
        uint32_t newDistance = matches[p->numPairs - 1];
        if ((p->longestMatchLength >= mainLen && newDistance < mainDist) ||
                (p->longestMatchLength == mainLen + 1 && !ChangePair(mainDist, newDistance)) ||
                (p->longestMatchLength > mainLen + 1) ||
                (p->longestMatchLength + 1 >= mainLen && mainLen >= 3 && ChangePair(newDistance, mainDist)))
            return 1;
    }

    data = Mf_GetPointerToCurrentPos(&p->matchFinderBase) - 1;
    for (i = 0; i < LZMA_NUM_REPS; i++)
    {
        uint32_t len, limit;
        const uint8_t *data2 = data - (p->reps[i] + 1);
        if (data[0] != data2[0] || data[1] != data2[1])
            continue;
        limit = mainLen - 1;
        for (len = 2; len < limit && data[len] == data2[len]; len++) ;
        if (len >= limit)
            return 1;
    }
    *backRes = mainDist + LZMA_NUM_REPS;
    MovePos(p, mainLen - 2);
    return mainLen;
}

static void LZe_full_flush(CLzmaEnc *p, uint32_t posState)
    {
    const uint32_t len = LZMA_MATCH_LEN_MIN;
    File_trailer trailer;
    RangeEnc_EncodeBit(&p->rc, &p->isMatch[p->state][posState], 1);
    RangeEnc_EncodeBit(&p->rc, &p->isRep[p->state], 0);
    p->state = kMatchNextStates[p->state];
    LenEnc_Encode2(&p->lenEnc, &p->rc, len - LZMA_MATCH_LEN_MIN, posState, !p->fastMode, p->ProbPrices);
    RcTree_Encode(&p->rc, p->posSlotEncoder[GetLenToPosState(len)], kNumPosSlotBits, (1 << kNumPosSlotBits) - 1);
    RangeEnc_EncodeDirectBits(&p->rc, (((uint32_t)1 << 30) - 1) >> kNumAlignBits, 30 - kNumAlignBits);
    RcTree_ReverseEncode(&p->rc, p->posAlignEncoder, kNumAlignBits, kAlignMask);
    RangeEnc_FlushData(&p->rc);
    RangeEnc_FlushStream(&p->rc);
    Ft_set_data_crc( trailer, p->matchFinderBase.crc ^ 0xFFFFFFFFU );
    Ft_set_data_size( trailer, p->nowPos64 );
    Ft_set_member_size( trailer, p->rc.processed + Fh_size + Ft_size );
    if( writeblock( p->rc.outfd, trailer, Ft_size ) != Ft_size )
        p->rc.res = SZ_ERROR_WRITE;
    if( verbosity >= 1 )
        {
        unsigned long long in_size = p->nowPos64;
        unsigned long long out_size = p->rc.processed + Fh_size + Ft_size;
        if( in_size == 0 || out_size == 0 )
            fputs( " no data compressed.\n", stderr );
        else
            fprintf( stderr, "%6.3f:1, %5.2f%% ratio, %5.2f%% saved, "
                                             "%llu in, %llu out.\n",
                             (double)in_size / out_size,
                             ( 100.0 * out_size ) / in_size,
                             100.0 - ( ( 100.0 * out_size ) / in_size ),
                             in_size, out_size );
        }
    }

static int CheckErrors(CLzmaEnc *p)
{
    if (p->result != SZ_OK)
        return p->result;
    if (p->rc.res != SZ_OK)
        p->result = SZ_ERROR_WRITE;
    if (p->matchFinderBase.result != SZ_OK)
        p->result = SZ_ERROR_READ;
    if (p->result != SZ_OK)
        p->finished = true;
    return p->result;
}

static int Flush(CLzmaEnc *p, uint32_t nowPos)
{
    /* ReleaseMFStream(); */
    p->finished = true;
    LZe_full_flush(p, nowPos & p->pbMask);
    return CheckErrors(p);
}

static void FillAlignPrices(CLzmaEnc *p)
{
    uint32_t i;
    for (i = 0; i < kAlignTableSize; i++)
        p->alignPrices[i] = RcTree_ReverseGetPrice(p->posAlignEncoder, kNumAlignBits, i, p->ProbPrices);
    p->alignPriceCount = 0;
}

static void FillDistancesPrices(CLzmaEnc *p)
{
    uint32_t tempPrices[kNumFullDistances];
    uint32_t i, lenToPosState;
    for (i = kStartPosModelIndex; i < kNumFullDistances; i++)
    {
        uint32_t posSlot = GetPosSlot1(i);
        uint32_t footerBits = ((posSlot >> 1) - 1);
        uint32_t base = ((2 | (posSlot & 1)) << footerBits);
        tempPrices[i] = RcTree_ReverseGetPrice(p->posEncoders + base - posSlot - 1, footerBits, i - base, p->ProbPrices);
    }

    for (lenToPosState = 0; lenToPosState < kNumLenToPosStates; lenToPosState++)
    {
        uint32_t posSlot;
        const int *encoder = p->posSlotEncoder[lenToPosState];
        uint32_t *posSlotPrices = p->posSlotPrices[lenToPosState];
        for (posSlot = 0; posSlot < p->distTableSize; posSlot++)
            posSlotPrices[posSlot] = RcTree_GetPrice(encoder, kNumPosSlotBits, posSlot, p->ProbPrices);
        for (posSlot = kEndPosModelIndex; posSlot < p->distTableSize; posSlot++)
            posSlotPrices[posSlot] += ((((posSlot >> 1) - 1) - kNumAlignBits) << kNumBitPriceShiftBits);

        {
            uint32_t *distancesPrices = p->distancesPrices[lenToPosState];
            uint32_t i;
            for (i = 0; i < kStartPosModelIndex; i++)
                distancesPrices[i] = posSlotPrices[i];
            for (; i < kNumFullDistances; i++)
                distancesPrices[i] = posSlotPrices[GetPosSlot1(i)] + tempPrices[i];
        }
    }
    p->matchPriceCount = 0;
}


static int LzmaEnc_CodeOneBlock(CLzmaEnc *p)
{
    uint32_t nowPos32, startPos32;

    if (p->finished)
        return p->result;
    if( CheckErrors(p) != 0 ) return p->result;

    nowPos32 = (uint32_t)p->nowPos64;
    startPos32 = nowPos32;

    if (p->nowPos64 == 0)
    {
        uint32_t numPairs;
        uint8_t curByte;
        if (Mf_GetNumAvailableBytes(&p->matchFinderBase) == 0)
            return Flush(p, nowPos32);
        ReadMatchDistances(p, &numPairs);
        RangeEnc_EncodeBit(&p->rc, &p->isMatch[p->state][0], 0);
        p->state = kLiteralNextStates[p->state];
        curByte = Mf_GetIndexByte(&p->matchFinderBase, 0 - p->additionalOffset);
        LitEnc_Encode(&p->rc, p->litProbs, curByte);
        p->additionalOffset--;
        nowPos32++;
    }

    if (Mf_GetNumAvailableBytes(&p->matchFinderBase) != 0)
    for (;;)
    {
        uint32_t pos, len, posState;

        if (p->fastMode)
            len = GetOptimumFast(p, &pos);
        else
            len = GetOptimum(p, nowPos32, &pos);

        #ifdef SHOW_STAT2
        printf("\n pos = %4X,   len = %d   pos = %d", nowPos32, len, pos);
        #endif

        posState = nowPos32 & p->pbMask;
        if (len == 1 && pos == (uint32_t)-1)
        {
            uint8_t curByte;
            int *probs;
            const uint8_t *data;

            RangeEnc_EncodeBit(&p->rc, &p->isMatch[p->state][posState], 0);
            data = Mf_GetPointerToCurrentPos(&p->matchFinderBase) - p->additionalOffset;
            curByte = *data;
            probs = LIT_PROBS(nowPos32, *(data - 1));
            if (IsCharState(p->state))
                LitEnc_Encode(&p->rc, probs, curByte);
            else
                LitEnc_EncodeMatched(&p->rc, probs, curByte, *(data - p->reps[0] - 1));
            p->state = kLiteralNextStates[p->state];
        }
        else
        {
            RangeEnc_EncodeBit(&p->rc, &p->isMatch[p->state][posState], 1);
            if (pos < LZMA_NUM_REPS)
            {
                RangeEnc_EncodeBit(&p->rc, &p->isRep[p->state], 1);
                if (pos == 0)
                {
                    RangeEnc_EncodeBit(&p->rc, &p->isRepG0[p->state], 0);
                    RangeEnc_EncodeBit(&p->rc, &p->isRep0Long[p->state][posState], ((len == 1) ? 0 : 1));
                }
                else
                {
                    uint32_t distance = p->reps[pos];
                    RangeEnc_EncodeBit(&p->rc, &p->isRepG0[p->state], 1);
                    if (pos == 1)
                        RangeEnc_EncodeBit(&p->rc, &p->isRepG1[p->state], 0);
                    else
                    {
                        RangeEnc_EncodeBit(&p->rc, &p->isRepG1[p->state], 1);
                        RangeEnc_EncodeBit(&p->rc, &p->isRepG2[p->state], pos - 2);
                        if (pos == 3)
                            p->reps[3] = p->reps[2];
                        p->reps[2] = p->reps[1];
                    }
                    p->reps[1] = p->reps[0];
                    p->reps[0] = distance;
                }
                if (len == 1)
                    p->state = kShortRepNextStates[p->state];
                else
                {
                    LenEnc_Encode2(&p->repLenEnc, &p->rc, len - LZMA_MATCH_LEN_MIN, posState, !p->fastMode, p->ProbPrices);
                    p->state = kRepNextStates[p->state];
                }
            }
            else
            {
                uint32_t posSlot;
                RangeEnc_EncodeBit(&p->rc, &p->isRep[p->state], 0);
                p->state = kMatchNextStates[p->state];
                LenEnc_Encode2(&p->lenEnc, &p->rc, len - LZMA_MATCH_LEN_MIN, posState, !p->fastMode, p->ProbPrices);
                pos -= LZMA_NUM_REPS;
                GetPosSlot(pos, posSlot);
                RcTree_Encode(&p->rc, p->posSlotEncoder[GetLenToPosState(len)], kNumPosSlotBits, posSlot);

                if (posSlot >= kStartPosModelIndex)
                {
                    uint32_t footerBits = ((posSlot >> 1) - 1);
                    uint32_t base = ((2 | (posSlot & 1)) << footerBits);
                    uint32_t posReduced = pos - base;

                    if (posSlot < kEndPosModelIndex)
                        RcTree_ReverseEncode(&p->rc, p->posEncoders + base - posSlot - 1, footerBits, posReduced);
                    else
                    {
                        RangeEnc_EncodeDirectBits(&p->rc, posReduced >> kNumAlignBits, footerBits - kNumAlignBits);
                        RcTree_ReverseEncode(&p->rc, p->posAlignEncoder, kNumAlignBits, posReduced & kAlignMask);
                        p->alignPriceCount++;
                    }
                }
                p->reps[3] = p->reps[2];
                p->reps[2] = p->reps[1];
                p->reps[1] = p->reps[0];
                p->reps[0] = pos;
                p->matchPriceCount++;
            }
        }
        p->additionalOffset -= len;
        nowPos32 += len;
        if (p->additionalOffset == 0)
        {
            uint32_t processed;
            if (!p->fastMode)
            {
                if (p->matchPriceCount >= (1 << 7))
                    FillDistancesPrices(p);
                if (p->alignPriceCount >= kAlignTableSize)
                    FillAlignPrices(p);
            }
            if (Mf_GetNumAvailableBytes(&p->matchFinderBase) == 0)
                break;
            processed = nowPos32 - startPos32;
            if (processed >= (1 << 15))
            {
                p->nowPos64 += nowPos32 - startPos32;
                return CheckErrors(p);
            }
        }
    }
    p->nowPos64 += nowPos32 - startPos32;
    return Flush(p, nowPos32);
}


CLzmaEncHandle LzmaEnc_Init( const int dict_size, const int match_len_limit,
                                                         const int infd, const int outfd )
    {
    int i;
    const uint32_t beforeSize = kNumOpts;
    CLzmaEnc * const p = (CLzmaEnc *)LZMA_MALLOC(sizeof(CLzmaEnc));
    if( !p ) return 0;

    p->nowPos64 = 0;
    p->dictSize = dict_size;
    p->numFastBytes = match_len_limit;
    p->lc = literal_context_bits;
    p->lp = 0;
    p->pb = pos_state_bits;
    p->optimumEndIndex = 0;
    p->optimumCurrentIndex = 0;
    p->additionalOffset = 0;
    p->state = 0;
    p->result = SZ_OK;
    p->fastMode = false;
    p->finished = false;

    if (!Mf_Init(&p->matchFinderBase, infd, 16 + ( match_len_limit / 2 ), p->dictSize, beforeSize, p->numFastBytes, LZMA_MATCH_LEN_MAX))
        { LZMA_FREE( p ); return 0; }
    Mf_CreateVTable(&p->matchFinderBase, &p->matchFinder);

    LzmaEnc_FastPosInit(p->g_FastPos);
    LzmaEnc_InitPriceTables(p->ProbPrices);
    for (i = 0; i < kDicLogSizeMaxCompress; i++)
        if (p->dictSize <= ((uint32_t)1 << i))
            break;
    p->distTableSize = i * 2;
    if( !RangeEnc_Init( &p->rc, outfd ) ) { LZMA_FREE( p ); return 0; }
    p->litProbs = (int *)LZMA_MALLOC((0x300 << (p->lc + p->lp)) * sizeof(int));
    if( !p->litProbs ) { LZMA_FREE( p ); return 0; }

    for (i = 0 ; i < LZMA_NUM_REPS; i++)
        p->reps[i] = 0;
    for (i = 0; i < kNumStates; i++)
    {
        int j;
        for (j = 0; j < LZMA_NUM_PB_STATES_MAX; j++)
        {
            p->isMatch[i][j] = kProbInitValue;
            p->isRep0Long[i][j] = kProbInitValue;
        }
        p->isRep[i] = kProbInitValue;
        p->isRepG0[i] = kProbInitValue;
        p->isRepG1[i] = kProbInitValue;
        p->isRepG2[i] = kProbInitValue;
    }
    {
        const int num = 0x300 << (p->lp + p->lc);
        for (i = 0; i < num; i++)
            p->litProbs[i] = kProbInitValue;
    }
    for (i = 0; i < kNumLenToPosStates; i++)
    {
        int *probs = p->posSlotEncoder[i];
        uint32_t j;
        for (j = 0; j < (1 << kNumPosSlotBits); j++)
            probs[j] = kProbInitValue;
    }
    for (i = 0; i < kNumFullDistances - kEndPosModelIndex; i++)
        p->posEncoders[i] = kProbInitValue;
    LenEnc_Init(&p->lenEnc.p);
    LenEnc_Init(&p->repLenEnc.p);
    for (i = 0; i < (1 << kNumAlignBits); i++)
        p->posAlignEncoder[i] = kProbInitValue;
    p->pbMask = (1 << p->pb) - 1;
    p->lpMask = (1 << p->lp) - 1;

    if (!p->fastMode) { FillDistancesPrices(p); FillAlignPrices(p); }
    p->lenEnc.tableSize =
    p->repLenEnc.tableSize =
            p->numFastBytes + 1 - LZMA_MATCH_LEN_MIN;
    LenPriceEnc_UpdateTables(&p->lenEnc, 1 << p->pb, p->ProbPrices);
    LenPriceEnc_UpdateTables(&p->repLenEnc, 1 << p->pb, p->ProbPrices);
    return p;
    }


void LzmaEnc_Free(CLzmaEncHandle pp)
{
    CLzmaEnc *p = (CLzmaEnc *)pp;
    Mf_Free(&p->matchFinderBase);
    LZMA_FREE(p->litProbs);
    p->litProbs = 0;
    RangeEnc_Free(&p->rc);
    LZMA_FREE(p);
}


int LzmaEnc_Encode(CLzmaEncHandle pp)
{
    int res = SZ_OK;
    CLzmaEnc *p = (CLzmaEnc *)pp;

    for (;;)
    {
        res = LzmaEnc_CodeOneBlock(p);
        if( res != SZ_OK || p->finished )
            break;
    }
    return res;
}

/* LzmaDec.h -- LZMA Decoder
2009-02-07 : Igor Pavlov : Public domain */

/* ---------- LZMA Properties ---------- */

#define LZMA_PROPS_SIZE 5

/* ---------- LZMA Decoder state ---------- */

/* LZMA_REQUIRED_INPUT_MAX = number of required input bytes for worst case.
     Num bits = log2((2^11 / 31) ^ 22) + 26 < 134 + 26 = 160; */

#define LZMA_REQUIRED_INPUT_MAX 20

typedef struct
{
    int *probs;
    uint8_t *dic;
    const uint8_t *buf;
    uint32_t range, code;
    uint32_t dicPos;
    uint32_t dicBufSize;
    uint32_t processedPos;
    uint32_t checkDicSize;
    unsigned lc, lp, pb;
    State state;
    uint32_t reps[4];
    unsigned remainLen;
    uint32_t numProbs;
    unsigned tempBufSize;
    bool needFlush;
    uint8_t tempBuf[LZMA_REQUIRED_INPUT_MAX];
} CLzmaDec;


/* There are two types of LZMA streams:
         0) Stream with end mark. That end mark adds about 6 bytes to compressed size.
         1) Stream without end mark. You must know exact uncompressed size to decompress such stream. */

typedef enum
{
    LZMA_FINISH_ANY,   /* finish at any point */
    LZMA_FINISH_END    /* block must be finished at the end */
} ELzmaFinishMode;

/* ELzmaFinishMode has meaning only if the decoding reaches output limit !!!

     You must use LZMA_FINISH_END, when you know that current output buffer
     covers last bytes of block. In other cases you must use LZMA_FINISH_ANY.

     If LZMA decoder sees end marker before reaching output limit, it returns SZ_OK,
     and output value of destLen will be less than output buffer size limit.
     You can check status result also.

     You can use multiple checks to test data integrity after full decompression:
         1) Check Result and "status" variable.
         2) Check that output(destLen) = uncompressedSize, if you know real uncompressedSize.
         3) Check that output(srcLen) = compressedSize, if you know real compressedSize.
                You must use correct finish mode in that case. */

typedef enum
{
    LZMA_STATUS_NOT_SPECIFIED,               /* use main error code instead */
    LZMA_STATUS_FINISHED_WITH_MARK,          /* stream was finished with end mark. */
    LZMA_STATUS_NOT_FINISHED,                /* stream was not finished */
    LZMA_STATUS_NEEDS_MORE_INPUT,            /* you must provide more input bytes */
    LZMA_STATUS_MAYBE_FINISHED_WITHOUT_MARK  /* there is probability that stream was finished without end mark */
} ELzmaStatus;

/* ELzmaStatus is used only as output value for function call */


static bool LzmaDec_Init(CLzmaDec *p, const uint8_t *raw_props);
static void LzmaDec_Free(CLzmaDec *p);


/* ---------- Buffer Interface ---------- */

/* It's zlib-like interface.

finishMode:
    It has meaning only if the decoding reaches output limit (*destLen).
    LZMA_FINISH_ANY - Decode just destLen bytes.
    LZMA_FINISH_END - Stream must be finished after (*destLen).
*/

static bool LzmaDec_DecodeToBuf( CLzmaDec *p, uint8_t *dest, uint32_t *destLen,
                                                    const uint8_t *src, uint32_t *srcLen,
                                                    ELzmaFinishMode finishMode, ELzmaStatus *status );

/* LzmaDec.c -- LZMA Decoder
2009-09-20 : Igor Pavlov : Public domain */

#define kNumTopBits 24
#define kTopValue ((uint32_t)1 << kNumTopBits)

#define kNumBitModelTotalBits 11
#define kBitModelTotal (1 << kNumBitModelTotalBits)
#define kNumMoveBits 5

#define RC_INIT_SIZE 5

#define NORMALIZE if (range < kTopValue) { range <<= 8; code = (code << 8) | (*buf++); }

#define IF_BIT_0(p) ttt = *(p); NORMALIZE; bound = (range >> kNumBitModelTotalBits) * ttt; if (code < bound)
#define UPDATE_0(p) range = bound; *(p) = (int)(ttt + ((kBitModelTotal - ttt) >> kNumMoveBits));
#define UPDATE_1(p) range -= bound; code -= bound; *(p) = (int)(ttt - (ttt >> kNumMoveBits));
#define GET_BIT2(p, i, A0, A1) IF_BIT_0(p) \
    { UPDATE_0(p); i = (i + i); A0; } else \
    { UPDATE_1(p); i = (i + i) + 1; A1; }
#define GET_BIT(p, i) GET_BIT2(p, i, ; , ;)

#define TREE_GET_BIT(probs, i) { GET_BIT((probs + i), i); }
#define TREE_DECODE(probs, limit, i) \
    { i = 1; do { TREE_GET_BIT(probs, i); } while (i < limit); i -= limit; }

/* #define _LZMA_SIZE_OPT */

#ifdef _LZMA_SIZE_OPT
#define TREE_6_DECODE(probs, i) TREE_DECODE(probs, (1 << 6), i)
#else
#define TREE_6_DECODE(probs, i) \
    { i = 1; \
    TREE_GET_BIT(probs, i); \
    TREE_GET_BIT(probs, i); \
    TREE_GET_BIT(probs, i); \
    TREE_GET_BIT(probs, i); \
    TREE_GET_BIT(probs, i); \
    TREE_GET_BIT(probs, i); \
    i -= 0x40; }
#endif

#define NORMALIZE_CHECK if (range < kTopValue) { if (buf >= bufLimit) return DUMMY_ERROR; range <<= 8; code = (code << 8) | (*buf++); }

#define IF_BIT_0_CHECK(p) ttt = *(p); NORMALIZE_CHECK; bound = (range >> kNumBitModelTotalBits) * ttt; if (code < bound)
#define UPDATE_0_CHECK range = bound;
#define UPDATE_1_CHECK range -= bound; code -= bound;
#define GET_BIT2_CHECK(p, i, A0, A1) IF_BIT_0_CHECK(p) \
    { UPDATE_0_CHECK; i = (i + i); A0; } else \
    { UPDATE_1_CHECK; i = (i + i) + 1; A1; }
#define GET_BIT_CHECK(p, i) GET_BIT2_CHECK(p, i, ; , ;)
#define TREE_DECODE_CHECK(probs, limit, i) \
    { i = 1; do { GET_BIT_CHECK(probs + i, i) } while (i < limit); i -= limit; }


#define kNumPosBitsMax 4
#define kNumPosStatesMax (1 << kNumPosBitsMax)

#define kLenNumLowBits 3
#define kLenNumLowSymbols (1 << kLenNumLowBits)
#define kLenNumMidBits 3
#define kLenNumMidSymbols (1 << kLenNumMidBits)
#define kLenNumHighBits 8
#define kLenNumHighSymbols (1 << kLenNumHighBits)

#define LenChoice 0
#define LenChoice2 (LenChoice + 1)
#define LenLow (LenChoice2 + 1)
#define LenMid (LenLow + (kNumPosStatesMax << kLenNumLowBits))
#define LenHigh (LenMid + (kNumPosStatesMax << kLenNumMidBits))
#define kNumLenProbs (LenHigh + kLenNumHighSymbols)


#define kNumStates 12
#define kNumLitStates 7

#define kStartPosModelIndex 4
#define kEndPosModelIndex 14
#define kNumFullDistances (1 << (kEndPosModelIndex >> 1))

#define kNumPosSlotBits 6
#define kNumLenToPosStates 4

#define kNumAlignBits 4
#define kAlignTableSize (1 << kNumAlignBits)

#define kMatchMinLen 2
#define kMatchSpecLenStart (kMatchMinLen + kLenNumLowSymbols + kLenNumMidSymbols + kLenNumHighSymbols)

#define IsMatch 0
#define IsRep (IsMatch + (kNumStates << kNumPosBitsMax))
#define IsRepG0 (IsRep + kNumStates)
#define IsRepG1 (IsRepG0 + kNumStates)
#define IsRepG2 (IsRepG1 + kNumStates)
#define IsRep0Long (IsRepG2 + kNumStates)
#define PosSlot (IsRep0Long + (kNumStates << kNumPosBitsMax))
#define SpecPos (PosSlot + (kNumLenToPosStates << kNumPosSlotBits))
#define Align (SpecPos + kNumFullDistances - kEndPosModelIndex)
#define LenCoder (Align + kAlignTableSize)
#define RepLenCoder (LenCoder + kNumLenProbs)
#define Literal (RepLenCoder + kNumLenProbs)

#define LZMA_BASE_SIZE 1846
#define LZMA_LIT_SIZE 768

#define LzmaProps_GetNumProbs(p) ((uint32_t)LZMA_BASE_SIZE + (LZMA_LIT_SIZE << ((p)->lc + (p)->lp)))

#if Literal != LZMA_BASE_SIZE
StopCompilingDueBUG
#endif


/* First LZMA-symbol is always decoded.
And it decodes new LZMA-symbols while (buf < bufLimit), but "buf" is without last normalization
Out:
    Result:
        true - OK
        false - Error
    p->remainLen:
        < kMatchSpecLenStart : normal remain
        = kMatchSpecLenStart : finished
        = kMatchSpecLenStart + 1 : Flush marker
        = kMatchSpecLenStart + 2 : State Init Marker
*/

static bool LzmaDec_DecodeReal(CLzmaDec *p, uint32_t limit, const uint8_t *bufLimit)
{
    int *probs = p->probs;

    State state = p->state;
    uint32_t rep0 = p->reps[0], rep1 = p->reps[1], rep2 = p->reps[2], rep3 = p->reps[3];
    unsigned pbMask = ((unsigned)1 << (p->pb)) - 1;
    unsigned lpMask = ((unsigned)1 << (p->lp)) - 1;
    const unsigned lc = p->lc;

    uint8_t *dic = p->dic;
    const uint32_t dicBufSize = p->dicBufSize;
    uint32_t dicPos = p->dicPos;

    uint32_t processedPos = p->processedPos;
    uint32_t checkDicSize = p->checkDicSize;
    unsigned len = 0;

    const uint8_t *buf = p->buf;
    uint32_t range = p->range;
    uint32_t code = p->code;

    do
    {
        int *prob;
        uint32_t bound;
        unsigned ttt;
        unsigned posState = processedPos & pbMask;

        prob = probs + IsMatch + (state << kNumPosBitsMax) + posState;
        IF_BIT_0(prob)
        {
            unsigned symbol;
            UPDATE_0(prob);
            prob = probs + Literal;
            if (checkDicSize != 0 || processedPos != 0)
                prob += (LZMA_LIT_SIZE * (((processedPos & lpMask) << lc) +
                (dic[(dicPos == 0 ? dicBufSize : dicPos) - 1] >> (8 - lc))));

            if (state < kNumLitStates)
            {
                state -= (state < 4) ? state : 3;
                symbol = 1;
                do { GET_BIT(prob + symbol, symbol) } while (symbol < 0x100);
            }
            else
            {
                unsigned matchByte = p->dic[(dicPos - rep0) + ((dicPos < rep0) ? dicBufSize : 0)];
                unsigned offs = 0x100;
                state -= (state < 10) ? 3 : 6;
                symbol = 1;
                do
                {
                    unsigned bit;
                    int *probLit;
                    matchByte <<= 1;
                    bit = (matchByte & offs);
                    probLit = prob + offs + bit + symbol;
                    GET_BIT2(probLit, symbol, offs &= ~bit, offs &= bit)
                }
                while (symbol < 0x100);
            }
            dic[dicPos++] = (uint8_t)symbol;
            processedPos++;
            continue;
        }
        else
        {
            UPDATE_1(prob);
            prob = probs + IsRep + state;
            IF_BIT_0(prob)
            {
                UPDATE_0(prob);
                state += kNumStates;
                prob = probs + LenCoder;
            }
            else
            {
                UPDATE_1(prob);
                if (checkDicSize == 0 && processedPos == 0)
                    return false;
                prob = probs + IsRepG0 + state;
                IF_BIT_0(prob)
                {
                    UPDATE_0(prob);
                    prob = probs + IsRep0Long + (state << kNumPosBitsMax) + posState;
                    IF_BIT_0(prob)
                    {
                        UPDATE_0(prob);
                        dic[dicPos] = dic[(dicPos - rep0) + ((dicPos < rep0) ? dicBufSize : 0)];
                        dicPos++;
                        processedPos++;
                        state = state < kNumLitStates ? 9 : 11;
                        continue;
                    }
                    UPDATE_1(prob);
                }
                else
                {
                    uint32_t distance;
                    UPDATE_1(prob);
                    prob = probs + IsRepG1 + state;
                    IF_BIT_0(prob)
                    {
                        UPDATE_0(prob);
                        distance = rep1;
                    }
                    else
                    {
                        UPDATE_1(prob);
                        prob = probs + IsRepG2 + state;
                        IF_BIT_0(prob)
                        {
                            UPDATE_0(prob);
                            distance = rep2;
                        }
                        else
                        {
                            UPDATE_1(prob);
                            distance = rep3;
                            rep3 = rep2;
                        }
                        rep2 = rep1;
                    }
                    rep1 = rep0;
                    rep0 = distance;
                }
                state = state < kNumLitStates ? 8 : 11;
                prob = probs + RepLenCoder;
            }
            {
                unsigned limit, offset;
                int *probLen = prob + LenChoice;
                IF_BIT_0(probLen)
                {
                    UPDATE_0(probLen);
                    probLen = prob + LenLow + (posState << kLenNumLowBits);
                    offset = 0;
                    limit = (1 << kLenNumLowBits);
                }
                else
                {
                    UPDATE_1(probLen);
                    probLen = prob + LenChoice2;
                    IF_BIT_0(probLen)
                    {
                        UPDATE_0(probLen);
                        probLen = prob + LenMid + (posState << kLenNumMidBits);
                        offset = kLenNumLowSymbols;
                        limit = (1 << kLenNumMidBits);
                    }
                    else
                    {
                        UPDATE_1(probLen);
                        probLen = prob + LenHigh;
                        offset = kLenNumLowSymbols + kLenNumMidSymbols;
                        limit = (1 << kLenNumHighBits);
                    }
                }
                TREE_DECODE(probLen, limit, len);
                len += offset;
            }

            if (state >= kNumStates)
            {
                uint32_t distance;
                prob = probs + PosSlot +
                        ((len < kNumLenToPosStates ? len : kNumLenToPosStates - 1) << kNumPosSlotBits);
                TREE_6_DECODE(prob, distance);
                if (distance >= kStartPosModelIndex)
                {
                    unsigned posSlot = (unsigned)distance;
                    int numDirectBits = (int)(((distance >> 1) - 1));
                    distance = (2 | (distance & 1));
                    if (posSlot < kEndPosModelIndex)
                    {
                        distance <<= numDirectBits;
                        prob = probs + SpecPos + distance - posSlot - 1;
                        {
                            uint32_t mask = 1;
                            unsigned i = 1;
                            do
                            {
                                GET_BIT2(prob + i, i, ; , distance |= mask);
                                mask <<= 1;
                            }
                            while (--numDirectBits != 0);
                        }
                    }
                    else
                    {
                        numDirectBits -= kNumAlignBits;
                        do
                        {
                            NORMALIZE
                            range >>= 1;

                            {
                                uint32_t t;
                                code -= range;
                                t = (0 - ((uint32_t)code >> 31)); /* (uint32_t)((int)code >> 31) */
                                distance = (distance << 1) + (t + 1);
                                code += range & t;
                            }
                            /*
                            distance <<= 1;
                            if (code >= range)
                            {
                                code -= range;
                                distance |= 1;
                            }
                            */
                        }
                        while (--numDirectBits != 0);
                        prob = probs + Align;
                        distance <<= kNumAlignBits;
                        {
                            unsigned i = 1;
                            GET_BIT2(prob + i, i, ; , distance |= 1);
                            GET_BIT2(prob + i, i, ; , distance |= 2);
                            GET_BIT2(prob + i, i, ; , distance |= 4);
                            GET_BIT2(prob + i, i, ; , distance |= 8);
                        }
                        if (distance == (uint32_t)0xFFFFFFFF)
                        {
                            len += kMatchSpecLenStart;
                            state -= kNumStates;
                            break;
                        }
                    }
                }
                rep3 = rep2;
                rep2 = rep1;
                rep1 = rep0;
                rep0 = distance + 1;
                if (checkDicSize == 0)
                {
                    if (distance >= processedPos)
                        return false;
                }
                else if (distance >= checkDicSize)
                    return false;
                state = (state < kNumStates + kNumLitStates) ? kNumLitStates : kNumLitStates + 3;
            }

            len += kMatchMinLen;

            if (limit == dicPos)
                return false;
            {
                uint32_t rem = limit - dicPos;
                unsigned curLen = ((rem < len) ? (unsigned)rem : len);
                uint32_t pos = (dicPos - rep0) + ((dicPos < rep0) ? dicBufSize : 0);

                processedPos += curLen;

                len -= curLen;
                if (pos + curLen <= dicBufSize)
                {
                    uint8_t *dest = dic + dicPos;
                    ptrdiff_t src = (ptrdiff_t)pos - (ptrdiff_t)dicPos;
                    const uint8_t *lim = dest + curLen;
                    dicPos += curLen;
                    do
                        *(dest) = (uint8_t)*(dest + src);
                    while (++dest != lim);
                }
                else
                {
                    do
                    {
                        dic[dicPos++] = dic[pos];
                        if (++pos == dicBufSize)
                            pos = 0;
                    }
                    while (--curLen != 0);
                }
            }
        }
    }
    while (dicPos < limit && buf < bufLimit);
    NORMALIZE;
    p->buf = buf;
    p->range = range;
    p->code = code;
    p->remainLen = len;
    p->dicPos = dicPos;
    p->processedPos = processedPos;
    p->reps[0] = rep0;
    p->reps[1] = rep1;
    p->reps[2] = rep2;
    p->reps[3] = rep3;
    p->state = state;

    return true;
}

static void LzmaDec_WriteRem(CLzmaDec *p, uint32_t limit)
{
    if (p->remainLen != 0 && p->remainLen < kMatchSpecLenStart)
    {
        uint8_t *dic = p->dic;
        uint32_t dicPos = p->dicPos;
        const uint32_t dicBufSize = p->dicBufSize;
        unsigned len = p->remainLen;
        uint32_t rep0 = p->reps[0];
        if (limit - dicPos < len)
            len = (unsigned)(limit - dicPos);

        if (p->checkDicSize == 0 && dicBufSize - p->processedPos <= len)
            p->checkDicSize = dicBufSize;

        p->processedPos += len;
        p->remainLen -= len;
        while (len-- != 0)
        {
            dic[dicPos] = dic[(dicPos - rep0) + ((dicPos < rep0) ? dicBufSize : 0)];
            dicPos++;
        }
        p->dicPos = dicPos;
    }
}

static int LzmaDec_DecodeReal2(CLzmaDec *p, uint32_t limit, const uint8_t *bufLimit)
{
    const uint32_t dicBufSize = p->dicBufSize;
    do
    {
        uint32_t limit2 = limit;
        if (p->checkDicSize == 0)
        {
            uint32_t rem = dicBufSize - p->processedPos;
            if (limit - p->dicPos > rem)
                limit2 = p->dicPos + rem;
        }
        if( !LzmaDec_DecodeReal(p, limit2, bufLimit) ) return false;
        if (p->processedPos >= dicBufSize)
            p->checkDicSize = dicBufSize;
        LzmaDec_WriteRem(p, limit);
    }
    while (p->dicPos < limit && p->buf < bufLimit && p->remainLen < kMatchSpecLenStart);

    if (p->remainLen > kMatchSpecLenStart)
    {
        p->remainLen = kMatchSpecLenStart;
    }
    return true;
}

typedef enum
{
    DUMMY_ERROR, /* unexpected end of input stream */
    DUMMY_LIT,
    DUMMY_MATCH,
    DUMMY_REP
} ELzmaDummy;

static ELzmaDummy LzmaDec_TryDummy(const CLzmaDec *p, const uint8_t *buf, uint32_t inSize)
{
    uint32_t range = p->range;
    uint32_t code = p->code;
    const uint8_t *bufLimit = buf + inSize;
    int *probs = p->probs;
    State state = p->state;
    ELzmaDummy res;

    {
        int *prob;
        uint32_t bound;
        unsigned ttt;
        unsigned posState = (p->processedPos) & ((1 << p->pb) - 1);

        prob = probs + IsMatch + (state << kNumPosBitsMax) + posState;
        IF_BIT_0_CHECK(prob)
        {
            UPDATE_0_CHECK

            /* if (bufLimit - buf >= 7) return DUMMY_LIT; */

            prob = probs + Literal;
            if (p->checkDicSize != 0 || p->processedPos != 0)
                prob += (LZMA_LIT_SIZE *
                    ((((p->processedPos) & ((1 << (p->lp)) - 1)) << p->lc) +
                    (p->dic[(p->dicPos == 0 ? p->dicBufSize : p->dicPos) - 1] >> (8 - p->lc))));

            if (state < kNumLitStates)
            {
                unsigned symbol = 1;
                do { GET_BIT_CHECK(prob + symbol, symbol) } while (symbol < 0x100);
            }
            else
            {
                unsigned matchByte = p->dic[p->dicPos - p->reps[0] +
                        ((p->dicPos < p->reps[0]) ? p->dicBufSize : 0)];
                unsigned offs = 0x100;
                unsigned symbol = 1;
                do
                {
                    unsigned bit;
                    int *probLit;
                    matchByte <<= 1;
                    bit = (matchByte & offs);
                    probLit = prob + offs + bit + symbol;
                    GET_BIT2_CHECK(probLit, symbol, offs &= ~bit, offs &= bit)
                }
                while (symbol < 0x100);
            }
            res = DUMMY_LIT;
        }
        else
        {
            unsigned len;
            UPDATE_1_CHECK;

            prob = probs + IsRep + state;
            IF_BIT_0_CHECK(prob)
            {
                UPDATE_0_CHECK;
                state = 0;
                prob = probs + LenCoder;
                res = DUMMY_MATCH;
            }
            else
            {
                UPDATE_1_CHECK;
                res = DUMMY_REP;
                prob = probs + IsRepG0 + state;
                IF_BIT_0_CHECK(prob)
                {
                    UPDATE_0_CHECK;
                    prob = probs + IsRep0Long + (state << kNumPosBitsMax) + posState;
                    IF_BIT_0_CHECK(prob)
                    {
                        UPDATE_0_CHECK;
                        NORMALIZE_CHECK;
                        return DUMMY_REP;
                    }
                    else
                    {
                        UPDATE_1_CHECK;
                    }
                }
                else
                {
                    UPDATE_1_CHECK;
                    prob = probs + IsRepG1 + state;
                    IF_BIT_0_CHECK(prob)
                    {
                        UPDATE_0_CHECK;
                    }
                    else
                    {
                        UPDATE_1_CHECK;
                        prob = probs + IsRepG2 + state;
                        IF_BIT_0_CHECK(prob)
                        {
                            UPDATE_0_CHECK;
                        }
                        else
                        {
                            UPDATE_1_CHECK;
                        }
                    }
                }
                state = kNumStates;
                prob = probs + RepLenCoder;
            }
            {
                unsigned limit, offset;
                int *probLen = prob + LenChoice;
                IF_BIT_0_CHECK(probLen)
                {
                    UPDATE_0_CHECK;
                    probLen = prob + LenLow + (posState << kLenNumLowBits);
                    offset = 0;
                    limit = 1 << kLenNumLowBits;
                }
                else
                {
                    UPDATE_1_CHECK;
                    probLen = prob + LenChoice2;
                    IF_BIT_0_CHECK(probLen)
                    {
                        UPDATE_0_CHECK;
                        probLen = prob + LenMid + (posState << kLenNumMidBits);
                        offset = kLenNumLowSymbols;
                        limit = 1 << kLenNumMidBits;
                    }
                    else
                    {
                        UPDATE_1_CHECK;
                        probLen = prob + LenHigh;
                        offset = kLenNumLowSymbols + kLenNumMidSymbols;
                        limit = 1 << kLenNumHighBits;
                    }
                }
                TREE_DECODE_CHECK(probLen, limit, len);
                len += offset;
            }

            if (state < 4)
            {
                unsigned posSlot;
                prob = probs + PosSlot +
                        ((len < kNumLenToPosStates ? len : kNumLenToPosStates - 1) <<
                        kNumPosSlotBits);
                TREE_DECODE_CHECK(prob, 1 << kNumPosSlotBits, posSlot);
                if (posSlot >= kStartPosModelIndex)
                {
                    int numDirectBits = ((posSlot >> 1) - 1);

                    /* if (bufLimit - buf >= 8) return DUMMY_MATCH; */

                    if (posSlot < kEndPosModelIndex)
                    {
                        prob = probs + SpecPos + ((2 | (posSlot & 1)) << numDirectBits) - posSlot - 1;
                    }
                    else
                    {
                        numDirectBits -= kNumAlignBits;
                        do
                        {
                            NORMALIZE_CHECK
                            range >>= 1;
                            code -= range & (((code - range) >> 31) - 1);
                            /* if (code >= range) code -= range; */
                        }
                        while (--numDirectBits != 0);
                        prob = probs + Align;
                        numDirectBits = kNumAlignBits;
                    }
                    {
                        unsigned i = 1;
                        do
                        {
                            GET_BIT_CHECK(prob + i, i);
                        }
                        while (--numDirectBits != 0);
                    }
                }
            }
        }
    }
    NORMALIZE_CHECK;
    return res;
}


static void LzmaDec_InitRc(CLzmaDec *p, const uint8_t *data)
{
    p->code = ((uint32_t)data[1] << 24) | ((uint32_t)data[2] << 16) | ((uint32_t)data[3] << 8) | ((uint32_t)data[4]);
    p->range = 0xFFFFFFFF;
    p->needFlush = false;
}


static bool LzmaDec_DecodeToDic(CLzmaDec *p, uint32_t dicLimit,
                                                                const uint8_t *src, uint32_t *srcLen,
                                                                ELzmaFinishMode finishMode, ELzmaStatus *status)
{
    uint32_t inSize = *srcLen;
    (*srcLen) = 0;
    LzmaDec_WriteRem(p, dicLimit);

    *status = LZMA_STATUS_NOT_SPECIFIED;

    while (p->remainLen != kMatchSpecLenStart)
    {
            int checkEndMarkNow;

            if( p->needFlush )
            {
                for (; inSize > 0 && p->tempBufSize < RC_INIT_SIZE; (*srcLen)++, inSize--)
                    p->tempBuf[p->tempBufSize++] = *src++;
                if (p->tempBufSize < RC_INIT_SIZE)
                {
                    *status = LZMA_STATUS_NEEDS_MORE_INPUT;
                    return true;
                }
                if (p->tempBuf[0] != 0)
                    return false;

                LzmaDec_InitRc(p, p->tempBuf);
                p->tempBufSize = 0;
            }

            checkEndMarkNow = 0;
            if (p->dicPos >= dicLimit)
            {
                if (p->remainLen == 0 && p->code == 0)
                {
                    *status = LZMA_STATUS_MAYBE_FINISHED_WITHOUT_MARK;
                    return true;
                }
                if (finishMode == LZMA_FINISH_ANY)
                {
                    *status = LZMA_STATUS_NOT_FINISHED;
                    return true;
                }
                if (p->remainLen != 0)
                {
                    *status = LZMA_STATUS_NOT_FINISHED;
                    return false;
                }
                checkEndMarkNow = 1;
            }

            if (p->tempBufSize == 0)
            {
                uint32_t processed;
                const uint8_t *bufLimit;
                if (inSize < LZMA_REQUIRED_INPUT_MAX || checkEndMarkNow)
                {
                    int dummyRes = LzmaDec_TryDummy(p, src, inSize);
                    if (dummyRes == DUMMY_ERROR)
                    {
                        memcpy(p->tempBuf, src, inSize);
                        p->tempBufSize = (unsigned)inSize;
                        (*srcLen) += inSize;
                        *status = LZMA_STATUS_NEEDS_MORE_INPUT;
                        return true;
                    }
                    if (checkEndMarkNow && dummyRes != DUMMY_MATCH)
                    {
                        *status = LZMA_STATUS_NOT_FINISHED;
                        return false;
                    }
                    bufLimit = src;
                }
                else
                    bufLimit = src + inSize - LZMA_REQUIRED_INPUT_MAX;
                p->buf = src;
                if( !LzmaDec_DecodeReal2(p, dicLimit, bufLimit) )
                    return false;
                processed = (uint32_t)(p->buf - src);
                (*srcLen) += processed;
                src += processed;
                inSize -= processed;
            }
            else
            {
                unsigned rem = p->tempBufSize, lookAhead = 0;
                while (rem < LZMA_REQUIRED_INPUT_MAX && lookAhead < inSize)
                    p->tempBuf[rem++] = src[lookAhead++];
                p->tempBufSize = rem;
                if (rem < LZMA_REQUIRED_INPUT_MAX || checkEndMarkNow)
                {
                    int dummyRes = LzmaDec_TryDummy(p, p->tempBuf, rem);
                    if (dummyRes == DUMMY_ERROR)
                    {
                        (*srcLen) += lookAhead;
                        *status = LZMA_STATUS_NEEDS_MORE_INPUT;
                        return true;
                    }
                    if (checkEndMarkNow && dummyRes != DUMMY_MATCH)
                    {
                        *status = LZMA_STATUS_NOT_FINISHED;
                        return false;
                    }
                }
                p->buf = p->tempBuf;
                if( !LzmaDec_DecodeReal2(p, dicLimit, p->buf) )
                    return false;
                lookAhead -= (rem - (unsigned)(p->buf - p->tempBuf));
                (*srcLen) += lookAhead;
                src += lookAhead;
                inSize -= lookAhead;
                p->tempBufSize = 0;
            }
    }
    if (p->code == 0)
        *status = LZMA_STATUS_FINISHED_WITH_MARK;
    return (p->code == 0);
}

static 
bool LzmaDec_DecodeToBuf( CLzmaDec *p, uint8_t *dest, uint32_t *destLen,
                                                    const uint8_t *src, uint32_t *srcLen,
                                                    ELzmaFinishMode finishMode, ELzmaStatus *status )
{
    uint32_t outSize = *destLen;
    uint32_t inSize = *srcLen;
    *srcLen = *destLen = 0;
    for (;;)
    {
        uint32_t inSizeCur = inSize, outSizeCur, dicPos;
        ELzmaFinishMode curFinishMode;
        bool res;
        if (p->dicPos == p->dicBufSize)
            p->dicPos = 0;
        dicPos = p->dicPos;
        if (outSize > p->dicBufSize - dicPos)
        {
            outSizeCur = p->dicBufSize;
            curFinishMode = LZMA_FINISH_ANY;
        }
        else
        {
            outSizeCur = dicPos + outSize;
            curFinishMode = finishMode;
        }

        res = LzmaDec_DecodeToDic(p, outSizeCur, src, &inSizeCur, curFinishMode, status);
        src += inSizeCur;
        inSize -= inSizeCur;
        *srcLen += inSizeCur;
        outSizeCur = p->dicPos - dicPos;
        memcpy(dest, p->dic + dicPos, outSizeCur);
        dest += outSizeCur;
        outSize -= outSizeCur;
        *destLen += outSizeCur;
        if( !res )
            return false;
        if (outSizeCur == 0 || outSize == 0)
            return true;
    }
}


static void LzmaDec_Free(CLzmaDec *p)
{
    LZMA_FREE( p->dic );
    LZMA_FREE( p->probs );
}


static bool LzmaDec_Init(CLzmaDec *p, const uint8_t *raw_props)
{
    uint32_t i;
    uint8_t d = raw_props[0];

    p->lc = d % 9;
    d /= 9;
    p->pb = d / 5;
    p->lp = d % 5;

    p->dicBufSize = raw_props[1] | ((uint32_t)raw_props[2] << 8) |
                                    ((uint32_t)raw_props[3] << 16) | ((uint32_t)raw_props[4] << 24);
    if (p->dicBufSize < min_dictionary_size) p->dicBufSize = min_dictionary_size;

    p->numProbs = LzmaProps_GetNumProbs(p);
    p->probs = (int *)LZMA_MALLOC(p->numProbs * sizeof(int));
    if( !p->probs ) return false;
    p->dic = (uint8_t *)LZMA_MALLOC(p->dicBufSize);
    if (p->dic == 0)
        {
        LZMA_FREE( p->probs );
        return false;
        }
    p->dicPos = 0;
    p->needFlush = true;
    p->remainLen = 0;
    p->tempBufSize = 0;
    p->processedPos = 0;
    p->checkDicSize = 0;
    for( i = 0; i < p->numProbs; ++i ) p->probs[i] = kBitModelTotal >> 1;
    p->reps[0] = p->reps[1] = p->reps[2] = p->reps[3] = 1;
    p->state = 0;
    return true;
}

// glue.c

static
#ifdef _MSC_VER
__declspec(thread)
#else
__thread
#endif
struct {
    uint8_t *begin, *seek, *end;
}
memfd[2];

/* Returns the number of bytes really read.
     If (returned value < size) and (errno == 0), means EOF was reached.
*/
static int readblock( const int fd, uint8_t * buf, int size ) {
    int avail = (memfd[fd].end - memfd[fd].seek);
    if( size > avail ) size = avail;
    memcpy(buf, memfd[fd].seek, size);
    memfd[fd].seek += size;
    errno = 0;
    return size;
}

/* Returns the number of bytes really written.
     If (returned value < size), it is always an error.
*/
static int writeblock( const int fd, const uint8_t *buf, int size ) {
    int avail = (memfd[fd].end - memfd[fd].seek);
    if( size > avail ) size = avail;
    memcpy(memfd[fd].seek, buf, size);
    memfd[fd].seek += size;
    errno = 0;
    return size;
}

// Customized compression modes. 
// Lower modes are optimized for low-mem devices. Uber modes A-B-C require *lots of RAM*.
static const struct lzma_options {
    int dictionary_size;      /* [4 KiB .. 512 MiB] */
    int match_len_limit;      /* [5 .. 273] */
} 
lzma_mappings[] = {
    // lowmem+fastest modes
    { 1 << 12,   5 },       // 0 - 39973598 lzma 39.97% c:13.635s d:2.909s
    { 1 << 16,   6 },       // 1 - 34979790 lzma 34.98% c:19.151s d:2.427s
    { 1 << 19,   7 },       // 2 - 32881806 lzma 32.88% c:25.592s d:1.907s
    { 1 << 20,   8 },       // 3 - 31908622 lzma 31.91% c:32.189s d:1.827s
    { 3 << 19,  10 },       // 4 - 30704458 lzma 30.70% c:40.736s d:1.747s
    { 1 << 21,  16 },       // 5 - 28807777 lzma 28.81% c:55.690s d:1.645s
    { 3 << 20,  20 },       // 6 - 28100304 lzma 28.10% c:63.734s d:1.614s
    { 1 << 22,  28 },       // 7 - 27594705 lzma 27.59% c:72.234s d:1.604s
    { 1 << 23,  36 },       // 8 - 27051139 lzma 27.05% c:79.418s d:1.586s
    { 1 << 24,  68 },       // 9 - 26702913 lzma 26.70% c:87.800s d:1.573s
    { 3 << 23, 132 },       // A - 26667550 lzma 26.67% c:89.020s d:1.581s
    { 1 << 25, 273 },       // B - 26656366 lzma 26.66% c:89.586s d:1.607s
    { 1 << 26, 273 },       // C - 26656366 lzma 26.66% c:90.004s d:1.586s
    // himem+slowest modes
};

unsigned lzma_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags /*[0..9]*/) {
    uint8_t level = (uint8_t)(flags > 9 ? 9 : flags < 0 ? 0 : flags);

    int i = 0; memfd[i].begin = memfd[i].seek = memfd[i].end = (uint8_t*)in; memfd[i].end += inlen;
    int o = 1; memfd[o].begin = memfd[o].seek = memfd[o].end = (uint8_t*)out;  memfd[o].end += outlen;

    writeblock(o, &level, 1); // write 1-byte header

    struct lzma_options encoder_options = lzma_mappings[level];
    CLzmaEncHandle handle = LzmaEnc_Init( encoder_options.dictionary_size, encoder_options.match_len_limit, i, o );
    int ok = SZ_OK == LzmaEnc_Encode(handle);
    LzmaEnc_Free(handle);

    return ok ? (int)(memfd[o].seek - memfd[o].begin) : 0;
}

unsigned lzma_decode(const void *in_, unsigned inlen, void *out, unsigned outlen) {
    const uint8_t *in = (const uint8_t*)in_;

    // parse 1-byte header
    uint8_t level = *in++; --inlen;

    //  -d{N}:  set dictionary size - [12, 30], default: 23 (8MB)
    //  -fb{N}: set number of fast bytes - [5, 273], default: 128
    //  -mc{N}: set number of cycles for match finder
    //  -lc{N}: set number of literal context bits - [0, 8], default: 3
    //  -lp{N}: set number of literal pos bits - [0, 4], default: 0
    //  -pb{N}: set number of pos bits - [0, 4], default: 2
    //  -mf{MF_ID}: set Match Finder: [bt2, bt3, bt4, hc4], default: bt4
    #pragma pack(push,1)
    struct { uint8_t d /*d=lc/pb/lp*/; uint32_t dsize; uint64_t rawsize; } props = {0};
    #pragma pack(pop)
    props.d = 0x5D;
    props.dsize = lzma_mappings[level].dictionary_size;

    CLzmaDec dec;
    ELzmaStatus status;
    LzmaDec_Init(&dec, &props.d);
    uint32_t srcLen = (uint32_t)inlen, destLen = (uint32_t)outlen;
    bool ok = LzmaDec_DecodeToBuf(&dec, (uint8_t*)out, &destLen, in, &srcLen, LZMA_FINISH_ANY, &status);
    LzmaDec_Free(&dec);

    return (unsigned)(ok ? destLen : 0);
}

unsigned lzma_bounds(unsigned inlen, unsigned flags) { 
    return (unsigned)(inlen * 1.1) + 16; // @todo: check src
}

#endif // LZMA_C

#ifdef LZMA_DEMO
#pragma once
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level = 1;
    char out[128];
    unsigned outlen = lzma_encode(longcopy, strlen(longcopy)+1, out, 128, level );
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    unsigned unpacked = lzma_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // LZMA_DEMO

#line 1 "amalgamated_lzp1.c" 
/***********

Direct port of the old lzp1.c code to a single file header.

This is not the best way to make fast compressors on modern hardware
and this is by no means a modern competitive compressor.

Also, zlib licensed is not strictly public domain, but pretty close terms :o)

-----------

Copyright (c) 2019, @r-lyeh
Copyright (c) 1998-2012, Charles Bloom

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

*******************/

unsigned lzp1_encode(const void* in, unsigned inlen, void* out, unsigned outlen, unsigned flags);
unsigned lzp1_decode(const void* in, unsigned inlen, void* out, unsigned outlen);
unsigned lzp1_bounds(unsigned inlen, unsigned flags);


#ifdef LZP1_C
#pragma once
#include <stdint.h>

#define LZP1_BOUNDS(sz)    ((sz)+((sz)/8)+256)
#define LZP1_EXCESS        256

#define LZP1_HASH_SIZE  (1<<16)
#define LZP1_HASH(x,y,z) ((x ^ (y << 7) ^ (z<<11)) & 0xFFFF)

static int lzp1_encode_(const uint8_t *raw,int rawLen,uint8_t * comp,int compLen)
{
    uint8_t const *table[LZP1_HASH_SIZE];
    for(int ix=0;ix<LZP1_HASH_SIZE;ix++) table[ix] = raw;

    uint8_t *cp,*controlp;
    const uint8_t *rp,*endrp,*mp;
    int ix,control,controlb,ml;
    uint8_t literal;

    /** do the LZP **/

    rp = raw;   endrp = raw + rawLen;
    cp = comp;

    // store excess
    *cp++ = rawLen & 255;

    // seed four
    *cp++ = *rp++; *cp++ = *rp++; *cp++ = *rp++; *cp++ = *rp++;

    control = 0; controlp = cp++; controlb = 8;

/** the control-byte entry macro **/

#define ENC_SHIFT_CONTROL(bit)  if ( 0 ) ; else { control += control + bit; if ( --controlb == 0 ) {    *controlp = (uint8_t)control;   controlp    = cp++; control     = 0;    controlb    = 8;    } }

    while(rp < endrp) {

        ix = LZP1_HASH(rp[-1],rp[-2],rp[-3]);
        mp = table[ix]; table[ix] = rp;

        if ( *mp != *rp ) {
            literal = *rp++;

            ix = LZP1_HASH(rp[-1],rp[-2],rp[-3]);
            mp = table[ix]; table[ix] = rp;

            if ( *mp != *rp ) {
                ENC_SHIFT_CONTROL(0);   //flag two literals : 0
                *cp++ = literal;
                *cp++ = *rp++;  // pass a literal
            } else {
                ENC_SHIFT_CONTROL(1);   //flag literal then a match : 10
                ENC_SHIFT_CONTROL(0);
                *cp++ = literal;
                goto encode_match;
            }
        } else {
            ENC_SHIFT_CONTROL(1);   //flag a match with no literals : 11
            ENC_SHIFT_CONTROL(1);

        encode_match:

            mp++; rp++;

            if ( *mp != *rp ) {
                ENC_SHIFT_CONTROL(0);
            } else {
                mp++; rp++;
                ENC_SHIFT_CONTROL(1);

                if ( *mp == *rp ) { mp++; rp++;
                    if ( *mp == *rp ) { mp++; rp++;
                        if ( *mp == *rp ) { mp++; rp++;
                            // flag more than 3
                            ENC_SHIFT_CONTROL(1);
                            ENC_SHIFT_CONTROL(1);
                            if ( *mp == *rp ) { mp++; rp++;
                                if ( *mp == *rp ) { mp++; rp++;
                                    if ( *mp == *rp ) { mp++; rp++;
                                        if ( *mp == *rp ) { mp++; rp++;
                                            if ( *mp == *rp ) { mp++; rp++;
                                                if ( *mp == *rp ) { mp++; rp++;
                                                    if ( *mp == *rp ) { mp++; rp++;
                                                            // flag 11 or more
                                                            ENC_SHIFT_CONTROL(1);
                                                            ENC_SHIFT_CONTROL(1);
                                                            ENC_SHIFT_CONTROL(1);

                                                            ml = 0;
                                                            while(rp < endrp  && *mp == *rp ) {
                                                                mp++; rp++; ml++;
                                                            }
                                                            while( ml >= 0xFF ) {
                                                                *cp++ = 0xFF; ml -= 0xFF;
                                                            }
                                                            *cp++ = (uint8_t)ml;

                                                    }   else {  // match 10
                                                        ENC_SHIFT_CONTROL(1);
                                                        ENC_SHIFT_CONTROL(1);
                                                        ENC_SHIFT_CONTROL(0);
                                                    }
                                                } else {    // match 9
                                                    ENC_SHIFT_CONTROL(1);
                                                    ENC_SHIFT_CONTROL(0);
                                                    ENC_SHIFT_CONTROL(1);
                                                }
                                            } else {    // match 8
                                                ENC_SHIFT_CONTROL(1);
                                                ENC_SHIFT_CONTROL(0);
                                                ENC_SHIFT_CONTROL(0);
                                            }
                                        } else {    // match 7
                                            ENC_SHIFT_CONTROL(0);
                                            ENC_SHIFT_CONTROL(1);
                                            ENC_SHIFT_CONTROL(1);
                                        }
                                    } else {    // match 6
                                        ENC_SHIFT_CONTROL(0);
                                        ENC_SHIFT_CONTROL(1);
                                        ENC_SHIFT_CONTROL(0);
                                    }
                                } else {    // match 5
                                    ENC_SHIFT_CONTROL(0);
                                    ENC_SHIFT_CONTROL(0);
                                    ENC_SHIFT_CONTROL(1);
                                }
                            }   else {  // match 4
                                ENC_SHIFT_CONTROL(0);
                                ENC_SHIFT_CONTROL(0);
                                ENC_SHIFT_CONTROL(0);
                            }
                        } else {    // match 3
                            ENC_SHIFT_CONTROL(1);
                            ENC_SHIFT_CONTROL(0);
                        }
                    } else {    // match 2
                        ENC_SHIFT_CONTROL(0);
                        ENC_SHIFT_CONTROL(1);
                    }
                } else {    //match 1
                    ENC_SHIFT_CONTROL(0);
                    ENC_SHIFT_CONTROL(0);
                }
            }
        }
    }

    //flush the control

    while( controlb > 0 ) {
        control += control;
        controlb--;
    }
    *controlp = (uint8_t)control;

    return (int)(cp - comp);
}

static int lzp1_decode_(const uint8_t * comp,int compLen,uint8_t * raw,int rawLen)
{
    uint8_t const *table[LZP1_HASH_SIZE];
    for(int ix=0;ix<LZP1_HASH_SIZE;ix++) table[ix] = raw;

    const uint8_t *cp,*mp,*endcp;
    uint8_t *rp,*endrp;
    int ix,control,controlb,ml;
    int bit;

    rp = raw;   endrp = raw + rawLen;
    cp = comp;  endcp = comp + compLen;

    uint8_t excess = *cp++; compLen--;

    *rp++ = *cp++; *rp++ = *cp++; *rp++ = *cp++; *rp++ = *cp++;

    control = *cp++;
    controlb = 8;

#define DEC_GET_CONTROL(getbit) if ( 0 ) ; else { getbit = control & 0x80;  control += control;if ( --controlb == 0 ) { control = *cp++;    controlb = 8; } }

    while(cp<endcp) {

        DEC_GET_CONTROL(bit);
        if ( ! bit ) {  // two literals
            table[ LZP1_HASH(rp[-1],rp[-2],rp[-3]) ] = rp;  *rp++ = *cp++;
            table[ LZP1_HASH(rp[-1],rp[-2],rp[-3]) ] = rp;  *rp++ = *cp++;

        } else {
            DEC_GET_CONTROL(bit);
            if ( ! bit ) {  //10 : literal then match
                table[ LZP1_HASH(rp[-1],rp[-2],rp[-3]) ] = rp;  *rp++ = *cp++;
            }   

            // match

            ix = LZP1_HASH(rp[-1],rp[-2],rp[-3]);   mp = table[ix]; table[ix] = rp;

            *rp++ = *mp++;

            // read 1 bit
            DEC_GET_CONTROL(bit);
            if ( bit ) {
                *rp++ = *mp++;
                // read 2 bits to get length
                DEC_GET_CONTROL(bit);
                if ( bit ) {
                    *rp++ = *mp++; *rp++ = *mp++;
                    DEC_GET_CONTROL(bit);
                    if ( bit ) {
                        *rp++ = *mp++;
                        //read 3 more bits

                        DEC_GET_CONTROL(bit);   if ( bit ) {
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;
                            DEC_GET_CONTROL(bit);   if ( bit ) {
                                DEC_GET_CONTROL(bit);   if ( bit ) {    // 111
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;

                                    do {
                                        int l;
                                        l = ml = *cp++;
                                        while(l--)
                                            *rp++ = *mp++;
                                    } while( ml == 0xFF );

                                } else {                                                        // 110
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;
                                }
                            } else {
                                DEC_GET_CONTROL(bit);   if ( bit ) {    // 101
                                    *rp++ = *mp++;
                                } else {                                                        // 100
                                }
                            }
                        } else {
                            DEC_GET_CONTROL(bit);   if ( bit ) {
                                DEC_GET_CONTROL(bit);   if ( bit ) {    // 011
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;
                                } else {                                                        // 010
                                    *rp++ = *mp++;
                                    *rp++ = *mp++;
                                }
                            } else {
                                DEC_GET_CONTROL(bit);   if ( bit ) {    // 001
                                    *rp++ = *mp++;
                                } else {                                                        // 000
                                }
                            }
                        }
                    }
                } else {
                    DEC_GET_CONTROL(bit);
                    if ( bit ) {
                        *rp++ = *mp++;
                    }
                }
            }
        }
    }

    return (((int)(rp - raw) >> 8) << 8) | excess;
}


unsigned lzp1_encode(const void* in, unsigned inlen, void* out, unsigned outlen, unsigned flags) {
    return (unsigned)lzp1_encode_((const uint8_t*)in, (int)inlen, (uint8_t*)out, (int)outlen);
}
unsigned lzp1_decode(const void* in, unsigned inlen, void* out, unsigned outlen) {
    return (unsigned)lzp1_decode_((const uint8_t*)in, (int)inlen, (uint8_t*)out, (int)outlen);
}
unsigned lzp1_bounds(unsigned inlen, unsigned flags) { 
    return (unsigned)LZP1_BOUNDS(inlen);
}

#endif // LZP1_C


#ifdef LZP1_DEMO
#pragma once
int main(int argc, char** argv) {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    char out[128];
    int outlen = lzp1_encode(longcopy, strlen(longcopy)+1, out, 128);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, outlen);

    char redo[128 + 256];
    int unpacked = lzp1_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", outlen, unpacked, redo);
}
#define main main__
#endif // LZP1_DEMO



#line 1 "amalgamated_lzrw3a.c" 
// Author : Ross Williams. Date : 15-Jul-1991. Release : 1.
// Modified by @r-lyeh.
//
// This file contains an implementation of the LZRW3-A data compression
// algorithm in the C programming language.
//    1 Algorithm is free of patent problems. The algorithm has not been
//      patented (nor will it be) and is of the LZ77 class which is fairly
//      clear of patents.
//    2 This implementation in C is in the public domain.

unsigned lzrw3a_encode(const void* in, unsigned inlen, void* out, unsigned outlen, unsigned flags);
unsigned lzrw3a_decode(const void* in, unsigned inlen, void* out, unsigned outlen);
unsigned lzrw3a_bounds(unsigned inlen, unsigned flags);


#ifdef LZRW3A_C
#pragma once
#include <string.h>
#include <stdint.h>

#define MEM_REQ ( HASH_TABLE_LENGTH*sizeof(uint8_t *) + 16 ) // 16 = ALIGNMENT_FUDGE

#define FLAG_BYTES 4
#define FLAG_COMPRESS 0
#define FLAG_COPY     1

#define ALIGN_UP(X) ((((uintptr_t)X)+3)&~3)

#define MAX_RAW_ITEM (18)
#define MAX_RAW_GROUP (16*MAX_RAW_ITEM)
#define MAX_CMP_GROUP (2+16*2)
#define HASH_TABLE_LENGTH (4096)
#define HASH_TABLE_DEPTH_BITS (3)
#define PARTITION_LENGTH_BITS (12-HASH_TABLE_DEPTH_BITS)
#define PARTITION_LENGTH      (1<<PARTITION_LENGTH_BITS)
#define HASH_TABLE_DEPTH      (1<<HASH_TABLE_DEPTH_BITS )
#define HASH_MASK             (PARTITION_LENGTH-1)
#define DEPTH_MASK            (HASH_TABLE_DEPTH-1)
#define START_STRING_18 ((uint8_t *) "123456789012345678")

#define HASH(PTR) ( \
         (((40543*(((*(PTR))<<8)^((*((PTR)+1))<<4)^(*((PTR)+2))))>>4) & HASH_MASK) \
    << HASH_TABLE_DEPTH_BITS \
)

#define UPDATE_P(P_BASE,NEWPTR) \
{(P_BASE)[cycle++]=(NEWPTR); cycle&=DEPTH_MASK;}

#define UPDATE_I(I_BASE,NEWPTR) \
{hash[(I_BASE)+cycle++]=(NEWPTR); cycle&=DEPTH_MASK;}

#define ANY_HASH_INDEX (0)

static void lzrw3a_compress(uint8_t* p_wrk_mem, uint8_t* p_src_first, uint32_t src_len, uint8_t* p_dst_first, size_t* p_dst_len)
{
    uint8_t* p_src = p_src_first;
    uint8_t* p_dst = p_dst_first;

    uint8_t* p_src_post = p_src_first + src_len;
    uint8_t* p_dst_post = p_dst_first + src_len;
    uint8_t* p_src_max1 = p_src_first + src_len - MAX_RAW_ITEM;
    uint8_t* p_src_max16 = p_src_first + src_len - MAX_RAW_ITEM * 16;

    #define TOPWORD 0xFFFF0000
    uint8_t* p_control;
    uint32_t control = TOPWORD;

    uint8_t** hash = (uint8_t**)ALIGN_UP(p_wrk_mem);

    uint8_t** p_h1 = 0;
    uint8_t** p_h2 = 0;

    unsigned cycle = 0;

    *p_dst++ = FLAG_COMPRESS;
    {unsigned i; for (i = 2; i <= FLAG_BYTES; i++) *p_dst++ = 0; }

    p_control = p_dst; p_dst += 2;

    {unsigned i; uint8_t** p_h = hash;
        #define ZH *p_h++=START_STRING_18
        for (i = 0; i < 256; i++)
        {
            ZH; ZH; ZH; ZH;
            ZH; ZH; ZH; ZH;
            ZH; ZH; ZH; ZH;
            ZH; ZH; ZH; ZH;
        }
    }

    while (1)
    {
        uint8_t* p_ziv;
        unsigned unroll;
        unsigned index;
        uint8_t** p_h0;
        register unsigned d;
        register unsigned bestlen;
        register unsigned bestpos;

        if (p_dst > p_dst_post)
            goto overrun;

        unroll = 16;
        if (p_src > p_src_max16)
        {
            unroll = 1;
            if (p_src > p_src_max1)
            {
                if (p_src == p_src_post)
                    break;
                else
                {
                    p_h0 = &hash[ANY_HASH_INDEX];
                    goto literal;
                }
            }
        }

        begin_unrolled_loop:
        p_ziv = p_src;

        index = HASH(p_src);
        p_h0 = &hash[index];
        bestlen = 0;
        bestpos = 0;
        for (d = 0; d < HASH_TABLE_DEPTH; d++)
        {
            register uint8_t* s = p_src;
            register uint8_t* p = p_h0[d];
            register unsigned len;
            if (s[bestlen] == p[bestlen])
            {
                #define PS *p++!=*s++
                PS || PS || PS || PS || PS || PS || PS || PS || PS ||
                    PS || PS || PS || PS || PS || PS || PS || PS || PS || s++;
                len = s - p_src - 1;
                if (len > bestlen)
                {
                    bestpos = d;
                    bestlen = len;
                }
            }
        }

        if (bestlen < 3)
        {
            literal: *p_dst++ = *p_src++; control &= 0xFFFEFFFF;
            if (p_h2 != 0)
            {
                UPDATE_P(p_h2, p_ziv - 2);
            }
            p_h2 = p_h1; p_h1 = p_h0;
        }
        else
        {
            index += bestpos;
            *p_dst++ = ((index & 0xF00) >> 4) | (bestlen - 3);
            *p_dst++ = index & 0xFF;
            p_src += bestlen;

            if (p_h1 != 0)
            {
                if (p_h2 != 0)
                {
                    UPDATE_P(p_h2, p_ziv - 2); p_h2 = 0;
                }
                UPDATE_P(p_h1, p_ziv - 1); p_h1 = 0;
            }
            UPDATE_P(p_h0, p_ziv);
        }
        control >>= 1;

        _unrolled_loop: 
        if (--unroll) goto begin_unrolled_loop;
        if ((control & TOPWORD) == 0)
        {
            *p_control++ = control & 0xFF;
            *p_control = (control >> 8) & 0xFF;

            p_control = p_dst; p_dst += 2;

            control = TOPWORD;
        }
    }
    while (control & TOPWORD) control >>= 1;
    *p_control++ = control & 0xFF;
    *p_control++ = (control >> 8) & 0xFF;

    if (p_control == p_dst) p_dst -= 2;

    *p_dst_len = p_dst - p_dst_first;
    return;
    overrun:
    *p_dst_first = FLAG_COPY;
    memcpy(p_dst_first + FLAG_BYTES, p_src_first, src_len);
    *p_dst_len = src_len + FLAG_BYTES;
}

static void lzrw3a_decompress(uint8_t* p_wrk_mem, uint8_t* p_src_first, uint32_t  src_len, uint8_t* p_dst_first, size_t* p_dst_len)
{
    register uint8_t* p_src = p_src_first + FLAG_BYTES;
    register uint8_t* p_dst = p_dst_first;

    uint8_t* p_src_post = p_src_first + src_len;
    uint8_t* p_src_max16 = p_src_first + src_len - (MAX_CMP_GROUP - 2);

    uint8_t** hash = (uint8_t**)ALIGN_UP(p_wrk_mem);    register uint32_t control = 1;
    register unsigned literals = 0;
    unsigned cycle = 0;
    if (*p_src_first == FLAG_COPY)
    {
        memcpy(p_dst_first, p_src_first + FLAG_BYTES, src_len - FLAG_BYTES);
        *p_dst_len = src_len - FLAG_BYTES;
        return;
    }
    {unsigned i; uint8_t** p_h = hash;
    #define ZJ *p_h++=START_STRING_18
        for (i = 0; i < 256; i++)
        {
            ZJ; ZJ; ZJ; ZJ;
            ZJ; ZJ; ZJ; ZJ;
            ZJ; ZJ; ZJ; ZJ;
            ZJ; ZJ; ZJ; ZJ;
        }
    }

    while (p_src != p_src_post)
    {
        register unsigned unroll;
        if (control == 1)
        {
            control = 0x10000 | *p_src++;
            control |= (*p_src++) << 8;
        }

        unroll = p_src <= p_src_max16 ? 16 : 1;

        while (unroll--)
        {
            if (control & 1)
            {
                register uint8_t* p;
                register unsigned lenmt;
                register uint8_t* p_ziv = p_dst;
                register unsigned index;

                lenmt = *p_src++;
                index = ((lenmt & 0xF0) << 4) | *p_src++;
                p = hash[index];
                lenmt &= 0xF;

                *p_dst++ = *p++;
                *p_dst++ = *p++;
                *p_dst++ = *p++;
                while (lenmt--)
                    *p_dst++ = *p++;
                if (literals > 0)
                {
                    register uint8_t* r = p_ziv - literals;;
                    UPDATE_I(HASH(r), r);
                    if (literals == 2)
                    {
                        r++; UPDATE_I(HASH(r), r);
                    }
                    literals = 0;
                }
                UPDATE_I(index & (~DEPTH_MASK), p_ziv);
            }
            else
            {
                *p_dst++ = *p_src++;

                if (++literals == 3)
                {
                    register uint8_t* p = p_dst - 3;
                    UPDATE_I(HASH(p), p); literals = 2;
                }
            }

            control >>= 1;
        }
    }

    *p_dst_len = p_dst - p_dst_first;
}

unsigned lzrw3a_encode(const void* in, unsigned inlen, void* out, unsigned outlen, unsigned flags) {
    uint8_t workmem[MEM_REQ];
    size_t outlen_ = outlen;
    lzrw3a_compress(workmem, (uint8_t*)in, inlen, (uint8_t*)out, &outlen_);
    return (unsigned)outlen_;
}
unsigned lzrw3a_decode(const void* in, unsigned inlen, void* out, unsigned outlen) {
    uint8_t workmem[MEM_REQ];
    size_t outlen_ = outlen;
    lzrw3a_decompress(workmem, (uint8_t*)in, inlen, (uint8_t*)out, &outlen_);
    return (unsigned)outlen_;
}
unsigned lzrw3a_bounds(unsigned inlen, unsigned flags) { 
    return (unsigned)(inlen * 1.1) + 16; // @todo: check src
}

#endif // LZRW3A_C

#ifdef LZRW3A_DEMO
#pragma once
#include <stdio.h>
int main() {
    const char* longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level = 1;
    char out[128];
    size_t outlen = lzrw3a_encode(longcopy, strlen(longcopy) + 1, out, 128, level);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy) + 1, (int)outlen);

    char redo[128];
    size_t unpacked = lzrw3a_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // LZRW3A_DEMO

#line 1 "amalgamated_lzss.c" 
/**************************************************************
    LZSS.C -- A Data Compression Program
***************************************************************
     4/ 6/1989   Haruhiko Okumura
    30/12/2019   @r-lyeh
    Use, distribute, and modify this program freely.
**************************************************************/

unsigned lzss_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags);
unsigned lzss_decode(const void *in, unsigned inlen, void *out, unsigned outlen);
unsigned lzss_bounds(unsigned bytes, unsigned flags);


#ifdef LZSS_C
#pragma once
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define N         4096    /* size of ring buffer */
#define F           18    /* upper limit for match_length */
#define THRESHOLD    2    /* encode string into position and length if match_length is greater than this */
#define NIL          N    /* index for root of binary search trees */

/* of longest match.  These are set by the InsertNode() procedure. */
static int match_position;
static int match_length;

static void InsertNode(unsigned char* text_buf, int* lson, int* rson, int* dad, int r)
    /* Inserts string of length F, text_buf[r..r+F-1], into one of the
       trees (text_buf[r]'th tree) and returns the longest-match position
       and length via the global variables match_position and match_length.
       If match_length = F, then removes the old node in favor of the new
       one, because the old one will be deleted sooner.
       Note r plays double role, as tree node and position in buffer. */
{
    int  i, p, cmp;
    unsigned char  *key;

    cmp = 1;  key = &text_buf[r];  p = N + 1 + key[0];
    rson[r] = lson[r] = NIL;  match_length = 0;
    for ( ; ; ) {
        if (cmp >= 0) {
            if (rson[p] != NIL) p = rson[p];
            else {  rson[p] = r;  dad[r] = p;  return;  }
        } else {
            if (lson[p] != NIL) p = lson[p];
            else {  lson[p] = r;  dad[r] = p;  return;  }
        }
        for (i = 1; i < F; i++)
            if ((cmp = key[i] - text_buf[p + i]) != 0)  break;
        if (i > match_length) {
            match_position = p;
            if ((match_length = i) >= F)  break;
        }
    }
    dad[r] = dad[p];  lson[r] = lson[p];  rson[r] = rson[p];
    dad[lson[p]] = r;  dad[rson[p]] = r;
    if (rson[dad[p]] == p) rson[dad[p]] = r;
    else                  lson[dad[p]] = r;
    dad[p] = NIL;  /* remove p */
}

static void DeleteNode(int* lson, int* rson, int* dad, int p)  /* deletes node p from tree */
{
    int  q;
    
    if (dad[p] == NIL) return;  /* not in tree */
    if (rson[p] == NIL) q = lson[p];
    else if (lson[p] == NIL) q = rson[p];
    else {
        q = lson[p];
        if (rson[q] != NIL) {
            do {  q = rson[q];  } while (rson[q] != NIL);
            rson[dad[q]] = lson[q];  dad[lson[q]] = dad[q];
            lson[q] = lson[p];  dad[lson[p]] = q;
        }
        rson[q] = rson[p];  dad[rson[p]] = q;
    }
    dad[q] = dad[p];
    if (rson[dad[p]] == p) rson[dad[p]] = q;  else lson[dad[p]] = q;
    dad[p] = NIL;
}

#define _get(c) \
    if (! ilen) {\
        c = -1; /*EOF;*/ \
        break;\
    }\
    c = *istr;\
    ++istr;\
    --ilen

#define _put(c) \
    *ostr = c;\
    ++ostr;\
    --olen

size_t LzssEncode(const char* istr, size_t ilen, char* ostr, size_t olen)
{
    int  i, c, len, r, s, last_match_length, code_buf_ptr;
    unsigned char  code_buf[17], mask;
    size_t codesize = 0;
    int lson[N + 1], rson[N + 257], dad[N + 1]; /* left & right children & parents -- These constitute binary search trees. */
    unsigned char text_buf[N + F - 1];          /* ring buffer of size N, with extra F-1 bytes to facilitate string comparison */

    match_position = 0;
    match_length = 0;

    if (ilen == 0) return 0;
    
    /* initialize trees */
    /* For i = 0 to N - 1, rson[i] and lson[i] will be the right and
    left children of node i.  These nodes need not be initialized.
    Also, dad[i] is the parent of node i.  These are initialized to
    NIL (= N), which stands for 'not used.'
    For i = 0 to 255, rson[N + i + 1] is the root of the tree
    for strings that begin with character i.  These are initialized
    to NIL.  Note there are 256 trees. */
    for (i = N + 1; i <= N + 256; i++) rson[i] = NIL;
    for (i = 0; i < N; i++) dad[i] = NIL;

    code_buf[0] = 0;  /* code_buf[1..16] saves eight units of code, and
        code_buf[0] works as eight flags, "1" representing that the unit
        is an unencoded letter (1 byte), "0" a position-and-length pair
        (2 bytes).  Thus, eight units require at most 16 bytes of code. */
    code_buf_ptr = mask = 1;
    s = 0;  r = N - F;
    for (i = s; i < r; i++) text_buf[i] = 0;  /* Clear the buffer with
        any character that will appear often. */
    for (len = 0; len < F && ilen; len++) {
        _get(c);
        text_buf[r + len] = c;
        /* Read F bytes into the last F bytes of the buffer */
    }
    for (i = 1; i <= F; i++) InsertNode(text_buf, lson, rson, dad, r - i);  /* Insert the F strings,
        each of which begins with one or more 'space' characters.  Note
        the order in which these strings are inserted.  This way,
        degenerate trees will be less likely to occur. */
    InsertNode(text_buf, lson, rson, dad, r);  /* Finally, insert the whole string just read. The global variables match_length and match_position are set. */
    do {
        if (match_length > len) match_length = len;  /* match_length may be spuriously long near the end of text. */
        if (match_length <= THRESHOLD) {
            match_length = 1;  /* Not long enough match.  Send one byte. */
            code_buf[0] |= mask;  /* 'send one byte' flag */
            code_buf[code_buf_ptr++] = text_buf[r];  /* Send uncoded. */
        } else {
            code_buf[code_buf_ptr++] = (unsigned char) match_position;
            code_buf[code_buf_ptr++] = (unsigned char)
                (((match_position >> 4) & 0xf0)
             | (match_length - (THRESHOLD + 1)));  /* Send position and
                    length pair. Note match_length > THRESHOLD. */
        }
        if ((mask <<= 1) == 0) {  /* Shift mask left one bit. */
            for (i = 0; i < code_buf_ptr; i++) {  /* Send at most 8 units of */
                _put(code_buf[i]);    /* code together */
            }
            codesize += code_buf_ptr;
            code_buf[0] = 0;  code_buf_ptr = mask = 1;
        }
        last_match_length = match_length;
        for (i = 0; i < last_match_length && ilen; i++) {
            _get(c);
            DeleteNode(lson, rson, dad, s);        /* Delete old strings and */
            text_buf[s] = c;    /* read new bytes */
            if (s < F - 1) text_buf[s + N] = c;  /* If the position is
                near the end of buffer, extend the buffer to make
                string comparison easier. */
            s = (s + 1) & (N - 1);  r = (r + 1) & (N - 1);
                /* Since this is a ring buffer, increment the position
                   modulo N. */
            InsertNode(text_buf, lson, rson, dad, r);    /* Register the string in text_buf[r..r+F-1] */
        }
        while (i++ < last_match_length) {    /* After the end of text, */
            DeleteNode(lson, rson, dad, s);                    /* no need to read, but */
            s = (s + 1) & (N - 1);  r = (r + 1) & (N - 1);
            if (--len) InsertNode(text_buf, lson, rson, dad, r);        /* buffer may not be empty. */
        }
    } while (len > 0);    /* until length of string to be processed is zero */
    if (code_buf_ptr > 1) {        /* Send remaining code. */
        for (i = 0; i < code_buf_ptr; i++) {
            _put(code_buf[i]);
        }
        codesize += code_buf_ptr;
    }

    return codesize;
}

#undef _put
#define _put(c) \
    *ostr++ = c;

size_t LzssDecode(const unsigned char* istr, size_t ilen, char *ostr, size_t olen)    /* Just the reverse of Encode(). */
{
    unsigned char text_buf[N + F - 1];          /* ring buffer of size N, with extra F-1 bytes to facilitate string comparison */
    int  i, j, k, r, c;
    unsigned int  flags;
    int limit = ilen;
    char *obak = ostr;
    
    for (i = 0; i < N - F; i++) text_buf[i] = 0;
    r = N - F;  flags = 0;
    for ( ; ; ) {
        if (((flags >>= 1) & 256) == 0) {
            _get(c);
            flags = c | 0xff00;        /* uses higher byte cleverly */
        }                            /* to count eight */
        if (flags & 1) {
            _get(c);
            _put(c);
            text_buf[r++] = c;  r &= (N - 1);
        } else {
            _get(i);
            _get(j);
            i |= ((j & 0xf0) << 4);  j = (j & 0x0f) + THRESHOLD;
            for (k = 0; k <= j; k++) {
                c = text_buf[(i + k) & (N - 1)];
                _put(c);
                text_buf[r++] = c;  r &= (N - 1);
            }
        }
    }
    return (size_t)(ostr - obak);
}

#undef _get
#undef _put
#undef N
#undef F
#undef THRESHOLD
#undef NIL

unsigned lzss_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags) {
    size_t rc = LzssEncode((const char*)in, (size_t)inlen, (char*)out, (size_t)outlen);
    return (unsigned)rc;
}
unsigned lzss_decode(const void *in, unsigned inlen, void *out, unsigned outlen) {
    size_t rc = LzssDecode((const unsigned char*)in, (size_t)inlen, (char*)out, (size_t)outlen);
    return (unsigned)rc;
}
unsigned lzss_bounds(unsigned bytes, unsigned flags) { 
    return (unsigned)(bytes * 1.5) + 16; // @todo: check src
}

#endif // LZSS_C


#ifdef LZSS_DEMO
#pragma once
#include <stdio.h>
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level=1;
    char out[128];
    size_t outlen = lzss_encode(longcopy, strlen(longcopy)+1, out, 128, level);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    size_t unpacked = lzss_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // LZSS_DEMO

#line 1 "amalgamated_ppp.c" 
// pred.c -- Original code by Dave Rand's rendition of the predictor algorithm.
// Updated by: Ian Donaldson, Carsten Bormann. Additional modifications by @r-lyeh.
//
// There are no license fees or costs associated with using the Predictor algorithm.
// Use the following code at your own risk.

unsigned ppp_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags);
unsigned ppp_decode(const void *in, unsigned inlen, void *out, unsigned outlen);
unsigned ppp_bounds(unsigned inlen, unsigned flags);


#ifdef PPP_C
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The following hash code is the heart of the algorithm:
 * It builds a sliding hash sum of the previous 3-and-a-bit
 * characters which will be used to index the guess table.
 * A better hash function would result in additional compression,
 * at the expense of time.
 */

// original. enwik8: 61.730.508 c:0.729s d:0.453s
//#define PPP_HASH_TYPE unsigned short
//#define PPP_HASH_TABLE (65536)
//#define PPP_HASH(x) Hash = (Hash << 4) ^ (x) //

// improved. enwik8: 58.769.363 c:0.772s d:0.490s 
#define PPP_HASH_TYPE unsigned int
#define PPP_HASH_TABLE (1<<18) // 256K
#define PPP_HASH(x) Hash = ((Hash * 160) ^ (x)) & (PPP_HASH_TABLE-1) // see: https://encode.su/threads/1025-PREDICTOR-algorithm

static int ppp_compress(const unsigned char *source, int slen, unsigned char *dest, int dlen) {
    PPP_HASH_TYPE Hash = 0;
    unsigned char GuessTable[PPP_HASH_TABLE] = {0};
    unsigned char *orgdest = dest;

    while (slen) {
        unsigned char *flagdest = dest++, flags = 0; /* All guess wrong initially */
        for (int bitmask=1, i=0; i < 8 && slen; i++, bitmask <<= 1) {
            if (GuessTable[Hash] != *source) {
                GuessTable[Hash] = *source;
                *dest++ = *source; /* Guess wrong, output char */
            } else {
                flags |= bitmask; /* Guess was right - don't output */
            }
            PPP_HASH(*source++);slen--;
        }
        *flagdest = flags;
    }
    return(dest - orgdest);
}

static int ppp_decompress(const unsigned char *source, int slen, unsigned char *dest, int dlen) {
    int final = 1;
    PPP_HASH_TYPE Hash = 0;
    unsigned char GuessTable[PPP_HASH_TABLE] = {0};
    unsigned char *orgdest = dest;
    while (slen >= 9) {
        unsigned char flags = *source++;
        for (int i=0, bitmask = 1; i < 8; i++, bitmask <<= 1) {
            if (!(flags & bitmask)) {
                GuessTable[Hash] = *source;     /* Guess wrong */
                *dest = *source++;          /* Read from source */
                slen--;
            } else {
                *dest = GuessTable[Hash];       /* Guess correct */
            }
            PPP_HASH(*dest++);
        }
        slen--;
    }
    while (final && slen > 0) {
        unsigned char flags = *source++;
        slen--;
        for (int i=0, bitmask = 1; i < 8; i++, bitmask <<= 1) {
            if (!(flags & bitmask)) {
                if (!slen)
                    break;  /* we seem to be really done -- cabo */
                GuessTable[Hash] = *source;     /* Guess wrong */
                *dest = *source++;          /* Read from source */
                slen--;
            } else {
                *dest = GuessTable[Hash];       /* Guess correct */
            }
            PPP_HASH(*dest++);
        }
    }
    return (dest - orgdest); // len
}

unsigned ppp_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags) {
    return (unsigned)ppp_compress((const unsigned char *)in, (int)inlen, (unsigned char *)out, (int)outlen);
}
unsigned ppp_decode(const void *in, unsigned inlen, void *out, unsigned outlen) {
    return (unsigned)ppp_decompress((const unsigned char *)in, (int)inlen, (unsigned char *)out, (int)outlen);
}
unsigned ppp_bounds(unsigned inlen, unsigned flags) { 
    return inlen/8*9+9;
}

#endif // PPP_C

#ifdef PPP_DEMO
#pragma once
#include <stdio.h>
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level = 0;
    char out[128];
    unsigned outlen = ppp_encode(longcopy, strlen(longcopy)+1, out, 128, level);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    unsigned unpacked = ppp_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", outlen, unpacked, redo);
}
#define main main__
#endif // PPP_DEMO

#line 1 "amalgamated_raw.c" 
// raw memcpy de/encoder
// - rlyeh, public domain

#ifndef RAW_H
#define RAW_H

unsigned raw_encode(const void *in, unsigned inlen, void *out, unsigned outcap, unsigned flags);
unsigned raw_decode(const void *in, unsigned inlen, void *out, unsigned outcap);
unsigned raw_bounds(unsigned bytes, unsigned flags);

#endif

#ifdef RAW_C
#pragma once
#include <string.h>

unsigned raw_encode(const void *in, unsigned inlen, void *out, unsigned outcap, unsigned flags) {
    return memcpy(out, in, inlen), inlen;
}
unsigned raw_decode(const void *in, unsigned inlen, void *out, unsigned outcap) {
    return memcpy(out, in, inlen), inlen;
}
unsigned raw_bounds(unsigned bytes, unsigned flags) {
    return bytes;
}

#endif

#line 1 "amalgamated_ulz.c" 
// ULZ.HPP - An ultra-fast LZ77 compressor
// Original C++ code written and placed in the public domain by Ilya Muravyov (UNLICENSED)
// Modified by r-lyeh (UNLICENSED)

unsigned ulz_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags); // [0..(6)..9]
unsigned ulz_decode(const void *in, unsigned inlen, void *out, unsigned outlen);
unsigned ulz_bounds(unsigned inlen, unsigned flags);


#ifdef ULZ_C
#pragma once
#include <stdlib.h>
#include <stdint.h>
#ifndef ULZ_REALLOC
#define ULZ_REALLOC realloc
#endif

enum { 
    ULZ_EXCESS=16,

    ULZ_WINDOW_BITS=17, // Hard-coded
    ULZ_WINDOW_SIZE=1<<ULZ_WINDOW_BITS,
    ULZ_WINDOW_MASK=ULZ_WINDOW_SIZE-1,

    ULZ_MIN_MATCH=4,

    ULZ_HASH_BITS=19,
    ULZ_HASH_SIZE=1<<ULZ_HASH_BITS,
    ULZ_NIL=-1,
};

typedef struct ULZ_WORKMEM {
    int HashTable[ULZ_HASH_SIZE];
    int Prev[ULZ_WINDOW_SIZE];
} ULZ_WORKMEM;

// Utils

static inline uint16_t UnalignedLoad16(const void* p) {
    return *(const uint16_t*)(p);
}

static inline uint32_t UnalignedLoad32(const void* p) {
    return *(const uint32_t*)(p);
}

static inline void UnalignedStore16(void* p, uint16_t x) {
    *(uint16_t*)(p)=x;
}

static inline void UnalignedCopy64(void* d, const void* s) {
    *(uint64_t*)(d)=*(const uint64_t*)(s);
}

static inline void WildCopy(uint8_t* d, const uint8_t* s, int n) {
    UnalignedCopy64(d, s);

    for (int i=8; i<n; i+=8)
        UnalignedCopy64(d+i, s+i);
}

static inline uint32_t Hash32(const void* p) {
    return (UnalignedLoad32(p)*0x9E3779B9)>>(32-ULZ_HASH_BITS);
}

static inline void EncodeMod(uint8_t** p, uint32_t x) {
    while (x>=128) {
        x-=128;
        *(*p)++=128+(x&127);
        x>>=7;
    }
    *(*p)++=x;
}

static inline uint32_t DecodeMod(const uint8_t** p) {
    uint32_t x=0;
    for (int i=0; i<=21; i+=7) {
        const uint32_t c=*(*p)++;
        x+=c<<i;
        if (c<128)
            break;
    }
    return x;
}

// LZ77

static int UlzCompressFast(const uint8_t* in, int inlen, uint8_t* out, int outlen) {
    ULZ_WORKMEM *u =(ULZ_WORKMEM*)ULZ_REALLOC(0, sizeof(ULZ_WORKMEM));

    for (int i=0; i<ULZ_HASH_SIZE; ++i)
        u->HashTable[i]=ULZ_NIL;

    uint8_t* op=out;
    int anchor=0;

    int p=0;
    while (p<inlen) {
        int best_len=0;
        int dist=0;

        const int max_match=inlen-p;
        if (max_match>=ULZ_MIN_MATCH) {
            const int limit=(p-ULZ_WINDOW_SIZE) > ULZ_NIL ? (p-ULZ_WINDOW_SIZE) : ULZ_NIL;

            const uint32_t h=Hash32(&in[p]);
            int s=u->HashTable[h];
            u->HashTable[h]=p;

            if (s>limit && UnalignedLoad32(&in[s])==UnalignedLoad32(&in[p])) {
                int len=ULZ_MIN_MATCH;
                while (len<max_match && in[s+len]==in[p+len])
                    ++len;

                best_len=len;
                dist=p-s;
            }
        }

        if (best_len==ULZ_MIN_MATCH && (p-anchor)>=(7+128))
            best_len=0;

        if (best_len>=ULZ_MIN_MATCH) {
            const int len=best_len-ULZ_MIN_MATCH;
            const int token=((dist>>12)&16)+(len < 15 ? len : 15);

            if (anchor!=p) {
                const int run=p-anchor;
                if (run>=7) {
                    *op++=(7<<5)+token;
                    EncodeMod(&op, run-7);
                }
                else
                    *op++=(run<<5)+token;

                WildCopy(op, &in[anchor], run);
                op+=run;
            }
            else
                *op++=token;

            if (len>=15)
                EncodeMod(&op, len-15);

            UnalignedStore16(op, dist);
            op+=2;

            anchor=p+best_len;
            ++p;
            u->HashTable[Hash32(&in[p])]=p++;
            u->HashTable[Hash32(&in[p])]=p++;
            u->HashTable[Hash32(&in[p])]=p++;
            p=anchor;
        }
        else
            ++p;
    }

    if (anchor!=p) {
        const int run=p-anchor;
        if (run>=7) {
            *op++=7<<5;
            EncodeMod(&op, run-7);
        }
        else
            *op++=run<<5;

        WildCopy(op, &in[anchor], run);
        op+=run;
    }

    ULZ_REALLOC(u, 0);
    return op-out;
}

static int UlzCompress(const uint8_t* in, int inlen, uint8_t* out, int outlen, int level) {
    if (level<1 || level>9)
        return 0;
    const int max_chain=(level<9)?1<<level:1<<13;

    ULZ_WORKMEM *u = (ULZ_WORKMEM*)ULZ_REALLOC(0, sizeof(ULZ_WORKMEM));
    for (int i=0; i<ULZ_HASH_SIZE; ++i)
        u->HashTable[i]=ULZ_NIL;

    uint8_t* op=out;
    int anchor=0;

    int p=0;
    while (p<inlen) {
        int best_len=0;
        int dist=0;

        const int max_match=inlen-p;
        if (max_match>=ULZ_MIN_MATCH) {
            const int limit=(p-ULZ_WINDOW_SIZE) > ULZ_NIL ? (p-ULZ_WINDOW_SIZE) : ULZ_NIL;
            int chainlen=max_chain;

            int s=u->HashTable[Hash32(&in[p])];
            while (s>limit) {
                if (in[s+best_len]==in[p+best_len]
                        && UnalignedLoad32(&in[s])==UnalignedLoad32(&in[p])) {
                    int len=ULZ_MIN_MATCH;
                    while (len<max_match && in[s+len]==in[p+len])
                        ++len;

                    if (len>best_len) {
                        best_len=len;
                        dist=p-s;

                        if (len==max_match)
                            break;
                    }
                }

                if (--chainlen==0)
                    break;

                s=u->Prev[s&ULZ_WINDOW_MASK];
            }
        }

        if (best_len==ULZ_MIN_MATCH && (p-anchor)>=(7+128))
            best_len=0;

        if (level>=5 && best_len>=ULZ_MIN_MATCH && best_len<max_match
                && (p-anchor)!=6) {
            const int x=p+1;
            const int target_len=best_len+1;

            const int limit=(x-ULZ_WINDOW_SIZE) > ULZ_NIL ? (x-ULZ_WINDOW_SIZE) : ULZ_NIL;
            int chainlen=max_chain;

            int s=u->HashTable[Hash32(&in[x])];
            while (s>limit) {
                if (in[s+best_len]==in[x+best_len]
                        && UnalignedLoad32(&in[s])==UnalignedLoad32(&in[x])) {
                    int len=ULZ_MIN_MATCH;
                    while (len<target_len && in[s+len]==in[x+len])
                        ++len;

                    if (len==target_len) {
                        best_len=0;
                        break;
                    }
                }

                if (--chainlen==0)
                    break;

                s=u->Prev[s&ULZ_WINDOW_MASK];
            }
        }

        if (best_len>=ULZ_MIN_MATCH) {
            const int len=best_len-ULZ_MIN_MATCH;
            const int token=((dist>>12)&16)+(len < 15 ? len : 15);

            if (anchor!=p) {
                const int run=p-anchor;
                if (run>=7) {
                    *op++=(7<<5)+token;
                    EncodeMod(&op, run-7);
                }
                else
                    *op++=(run<<5)+token;

                WildCopy(op, &in[anchor], run);
                op+=run;
            }
            else
                *op++=token;

            if (len>=15)
                EncodeMod(&op, len-15);

            UnalignedStore16(op, dist);
            op+=2;

            while (best_len--!=0) {
                const uint32_t h=Hash32(&in[p]);
                u->Prev[p&ULZ_WINDOW_MASK]=u->HashTable[h];
                u->HashTable[h]=p++;
            }
            anchor=p;
        }
        else {
            const uint32_t h=Hash32(&in[p]);
            u->Prev[p&ULZ_WINDOW_MASK]=u->HashTable[h];
            u->HashTable[h]=p++;
        }
    }

    if (anchor!=p) {
        const int run=p-anchor;
        if (run>=7) {
            *op++=7<<5;
            EncodeMod(&op, run-7);
        }
        else
            *op++=run<<5;

        WildCopy(op, &in[anchor], run);
        op+=run;
    }

    ULZ_REALLOC(u, 0);
    return op-out;
}

static int UlzDecompress(const uint8_t* in, int inlen, uint8_t* out, int outlen) {
    uint8_t* op=out;
    const uint8_t* ip=in;
    const uint8_t* ip_end=ip+inlen;
    const uint8_t* op_end=op+outlen;

    while (ip<ip_end) {
        const int token=*ip++;

        if (token>=32) {
            int run=token>>5;
            if (run==7)
                run+=DecodeMod(&ip);
            if ((op_end-op)<run || (ip_end-ip)<run) // Overrun check
                return 0;

            WildCopy(op, ip, run);
            op+=run;
            ip+=run;
            if (ip>=ip_end)
                break;
        }

        int len=(token&15)+ULZ_MIN_MATCH;
        if (len==(15+ULZ_MIN_MATCH))
            len+=DecodeMod(&ip);
        if ((op_end-op)<len) // Overrun check
            return 0;

        const int dist=((token&16)<<12)+UnalignedLoad16(ip);
        ip+=2;
        uint8_t* cp=op-dist;
        if ((op-out)<dist) // Range check
            return 0;

        if (dist>=8) {
            WildCopy(op, cp, len);
            op+=len;
        }
        else {
            *op++=*cp++;
            *op++=*cp++;
            *op++=*cp++;
            *op++=*cp++;
            while (len--!=4)
                *op++=*cp++;
        }
    }

    return (ip==ip_end)?op-out:0;
}

unsigned ulz_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags) {
    int level = flags > 9 ? 9 : flags < 0 ? 0 : flags; // [0..(6)..9]
    int rc = level ? UlzCompress((uint8_t *)in, (int)inlen, (uint8_t *)out, (int)outlen, level)
        : UlzCompressFast((uint8_t *)in, (int)inlen, (uint8_t *)out, (int)outlen);
    return (unsigned)rc;
}
unsigned ulz_decode(const void *in, unsigned inlen, void *out, unsigned outlen) {
    return (unsigned)UlzDecompress((uint8_t *)in, (int)inlen, (uint8_t *)out, (int)outlen);
}
unsigned ulz_bounds(unsigned inlen, unsigned flags) { 
    return inlen + inlen/255 + 16;
}

#endif // ULZ_C


#ifdef ULZ_DEMO
#pragma once
#include <stdio.h>
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level=1;
    char out[128];
    size_t outlen = ulz_encode(longcopy, strlen(longcopy)+1, out, 128, level);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    size_t unpacked = ulz_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // ULZ_DEMO

#line 1 "amalgamated_pack.c" 
#ifdef STDPACK_C
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
    return 1 + list[(compressor >> 4) % NUM_COMPRESSORS].bounds(inlen, compressor & 0x0F);
}

// ---
// file options

static uint8_t STDPACK_FILE_BLOCK_SIZE =  23; // 2<<(BS+12) = { 8K..256M }
static uint8_t STDPACK_FILE_BLOCK_EXCESS = 0; // 16<<BE = 16, 256, 4K, 64K (16 for ulz, 256 for lpz1)

// xx yy zzzz   : 8 bits
// xx           : reserved (default = 0x11)
//    yy        : block excess [00..03] = 16<<X     = { 16, 256, 4K, 64K }
//       zzzz   : block size   [00..15] = 2<<(X+13) = { 8K..256M }

unsigned file_encode(FILE* in, FILE* out, FILE *logfile, unsigned cnum, unsigned *clist) { // multi encoder
#if 0
    // uint8_t MAGIC = 0x11 << 6 | ((STDPACK_FILE_BLOCK_EXCESS&3) << 4) | ((STDPACK_FILE_BLOCK_SIZE-12)&15);
    // EXCESS = 16ull << ((MAGIC >> 4) & 3);
    // BLSIZE =  1ull << ((MAGIC & 15) + 13);
#else
    if( fwrite(&STDPACK_FILE_BLOCK_SIZE, 1,1, out) < 1) return 0;
    if( fwrite(&STDPACK_FILE_BLOCK_EXCESS, 1,1, out) < 1) return 0;
    uint64_t BS_BYTES = 1ull << STDPACK_FILE_BLOCK_SIZE;
    uint64_t BE_BYTES = 1ull << STDPACK_FILE_BLOCK_EXCESS;
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

#endif // STDPACK_C

