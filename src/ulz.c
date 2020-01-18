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
