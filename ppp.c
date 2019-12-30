// pred.c -- Original code by Dave Rand's rendition of the predictor algorithm.
// Updated by: Ian Donaldson, Carsten Bormann. Additional modifications by @r-lyeh.
//
// There are no license fees or costs associated with using the Predictor algorithm.
// Use the following code at your own risk.

unsigned ppp_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags);
unsigned ppp_decode(const void *in, unsigned inlen, void *out, unsigned outlen);


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

// original
//#define PPP_HASH_TYPE unsigned short
//#define PPP_HASH_TABLE (65536)
//#define PPP_HASH(x) Hash = (Hash << 4) ^ (x) // 61.730.508 0.729s 0.453s

#define PPP_HASH_TYPE unsigned int
#define PPP_HASH_TABLE (1<<18) // 256K
#define PPP_HASH(x) Hash = ((Hash * 160) ^ (x)) & (PPP_HASH_TABLE-1) // 58.769.363 0.772s 0.490s // see: https://encode.su/threads/1025-PREDICTOR-algorithm

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
