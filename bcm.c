// BCM 1.40 - A BWT-based file compressor
// Written and placed in the public domain by Ilya Muravyov (UNLICENSE)
// Additional code by @r-lyeh (UNLICENSE)
//
// Notes:
// - BCM decoder has no dependencies.
// - BCM encoder requires libdivsufsort, which is MIT licensed.
// - #define BCM_NO_ENCODER if you want to exclude libdivsufsort from linkage.

unsigned bcm_bounds(unsigned inlen, unsigned flags);
unsigned bcm_encode(const void *in, unsigned inlen, const void *out, unsigned outlen, unsigned flags/*[0..(4)..9]*/);
unsigned bcm_decode(const void *in, unsigned inlen, const void *out, unsigned outlen);

// ---

#ifdef BCM_C
#pragma once
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef BCM_NO_ENCODER
int bcm_divbwt(const unsigned char *T, unsigned char *U, int *A, int n) { return -1; }
#else
#include "bcm_bwt.c" // libdivsufsort-lite (MIT licensed)
#endif

#if INTPTR_MAX >= INT64_MAX
#define BCM_64BITS 1
#else
#define BCM_64BITS 0
#endif

#ifndef BCM_REALLOC
#define BCM_REALLOC realloc
#endif

#if defined _MSC_VER && !defined __thread
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
bcmCounter##RATE##Ctor(bcmCounter##RATE *c) { c->p=1<<15; /* 0.5 */ } \
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
