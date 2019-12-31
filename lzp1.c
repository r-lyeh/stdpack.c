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
                                                            while( *mp == *rp ) {
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


