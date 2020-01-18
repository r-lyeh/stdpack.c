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

inline uint32_t get_hash(int p) {
	return (((*(uint32_t*)(&buf[p]))&0xffffff) *2654435769UL)&~BALZ_BUF_MASK; // Little-endian+unaligned
}

inline int get_pts(int len, int x) {
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
