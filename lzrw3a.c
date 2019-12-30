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

	uint8_t** hash = (uint8_t**)ALIGN_UP(p_wrk_mem);	register uint32_t control = 1;
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
	char workmem[MEM_REQ];
	size_t outlen_ = outlen;
	lzrw3a_compress(workmem, (void*)in, inlen, out, &outlen_);
	return (unsigned)outlen_;
}

unsigned lzrw3a_decode(const void* in, unsigned inlen, void* out, unsigned outlen) {
	char workmem[MEM_REQ];
	size_t outlen_ = outlen;
	lzrw3a_decompress(workmem, (void*)in, inlen, out, &outlen_);
	return (unsigned)outlen_;
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
