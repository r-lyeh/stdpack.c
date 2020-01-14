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
