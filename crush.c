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
