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
