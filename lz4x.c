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
  if(level >= 15) return lz4x_compress_optimal(in, inlen, out, outlen);
  return (unsigned)lz4x_compress((const uint8_t*)in, (size_t)inlen, (uint8_t*)out, (size_t)outlen, level);
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
