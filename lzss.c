/**************************************************************
    LZSS.C -- A Data Compression Program
    (tab = 4 spaces)
***************************************************************
    4/6/1989 Haruhiko Okumura
    Use, distribute, and modify this program freely.
    Please send me your improved versions.
        PC-VAN        SCIENCE
        NIFTY-Serve    PAF01022
        CompuServe    74050,1022
**************************************************************/

unsigned lzss_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags);
unsigned lzss_decode(const void *in, unsigned inlen, void *out, unsigned outlen);


#ifdef LZSS_C
#pragma once
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define N         4096    /* size of ring buffer */
#define F           18    /* upper limit for match_length */
#define THRESHOLD    2    /* encode string into position and length if match_length is greater than this */
#define NIL          N    /* index for root of binary search trees */

/* of longest match.  These are set by the InsertNode() procedure. */
static int match_position;
static int match_length;

static void InsertNode(unsigned char* text_buf, int* lson, int* rson, int* dad, int r)
    /* Inserts string of length F, text_buf[r..r+F-1], into one of the
       trees (text_buf[r]'th tree) and returns the longest-match position
       and length via the global variables match_position and match_length.
       If match_length = F, then removes the old node in favor of the new
       one, because the old one will be deleted sooner.
       Note r plays double role, as tree node and position in buffer. */
{
    int  i, p, cmp;
    unsigned char  *key;

    cmp = 1;  key = &text_buf[r];  p = N + 1 + key[0];
    rson[r] = lson[r] = NIL;  match_length = 0;
    for ( ; ; ) {
        if (cmp >= 0) {
            if (rson[p] != NIL) p = rson[p];
            else {  rson[p] = r;  dad[r] = p;  return;  }
        } else {
            if (lson[p] != NIL) p = lson[p];
            else {  lson[p] = r;  dad[r] = p;  return;  }
        }
        for (i = 1; i < F; i++)
            if ((cmp = key[i] - text_buf[p + i]) != 0)  break;
        if (i > match_length) {
            match_position = p;
            if ((match_length = i) >= F)  break;
        }
    }
    dad[r] = dad[p];  lson[r] = lson[p];  rson[r] = rson[p];
    dad[lson[p]] = r;  dad[rson[p]] = r;
    if (rson[dad[p]] == p) rson[dad[p]] = r;
    else                  lson[dad[p]] = r;
    dad[p] = NIL;  /* remove p */
}

static void DeleteNode(int* lson, int* rson, int* dad, int p)  /* deletes node p from tree */
{
    int  q;
    
    if (dad[p] == NIL) return;  /* not in tree */
    if (rson[p] == NIL) q = lson[p];
    else if (lson[p] == NIL) q = rson[p];
    else {
        q = lson[p];
        if (rson[q] != NIL) {
            do {  q = rson[q];  } while (rson[q] != NIL);
            rson[dad[q]] = lson[q];  dad[lson[q]] = dad[q];
            lson[q] = lson[p];  dad[lson[p]] = q;
        }
        rson[q] = rson[p];  dad[rson[p]] = q;
    }
    dad[q] = dad[p];
    if (rson[dad[p]] == p) rson[dad[p]] = q;  else lson[dad[p]] = q;
    dad[p] = NIL;
}

#define _get(c) \
    if (! ilen) {\
        c = -1; /*EOF;*/ \
        break;\
    }\
    c = *istr;\
    ++istr;\
    --ilen

#define _put(c) \
    *ostr = c;\
    ++ostr;\
    --olen

size_t LzssEncode(const char* istr, size_t ilen, char* ostr, size_t olen)
{
    int  i, c, len, r, s, last_match_length, code_buf_ptr;
    unsigned char  code_buf[17], mask;
    size_t codesize = 0;
    int lson[N + 1], rson[N + 257], dad[N + 1]; /* left & right children & parents -- These constitute binary search trees. */
    unsigned char text_buf[N + F - 1];          /* ring buffer of size N, with extra F-1 bytes to facilitate string comparison */

    match_position = 0;
    match_length = 0;

    if (ilen == 0) return 0;
    
    /* initialize trees */
    /* For i = 0 to N - 1, rson[i] and lson[i] will be the right and
    left children of node i.  These nodes need not be initialized.
    Also, dad[i] is the parent of node i.  These are initialized to
    NIL (= N), which stands for 'not used.'
    For i = 0 to 255, rson[N + i + 1] is the root of the tree
    for strings that begin with character i.  These are initialized
    to NIL.  Note there are 256 trees. */
    for (i = N + 1; i <= N + 256; i++) rson[i] = NIL;
    for (i = 0; i < N; i++) dad[i] = NIL;

    code_buf[0] = 0;  /* code_buf[1..16] saves eight units of code, and
        code_buf[0] works as eight flags, "1" representing that the unit
        is an unencoded letter (1 byte), "0" a position-and-length pair
        (2 bytes).  Thus, eight units require at most 16 bytes of code. */
    code_buf_ptr = mask = 1;
    s = 0;  r = N - F;
    for (i = s; i < r; i++) text_buf[i] = 0;  /* Clear the buffer with
        any character that will appear often. */
    for (len = 0; len < F && ilen; len++) {
        _get(c);
        text_buf[r + len] = c;
        /* Read F bytes into the last F bytes of the buffer */
    }
    for (i = 1; i <= F; i++) InsertNode(text_buf, lson, rson, dad, r - i);  /* Insert the F strings,
        each of which begins with one or more 'space' characters.  Note
        the order in which these strings are inserted.  This way,
        degenerate trees will be less likely to occur. */
    InsertNode(text_buf, lson, rson, dad, r);  /* Finally, insert the whole string just read. The global variables match_length and match_position are set. */
    do {
        if (match_length > len) match_length = len;  /* match_length may be spuriously long near the end of text. */
        if (match_length <= THRESHOLD) {
            match_length = 1;  /* Not long enough match.  Send one byte. */
            code_buf[0] |= mask;  /* 'send one byte' flag */
            code_buf[code_buf_ptr++] = text_buf[r];  /* Send uncoded. */
        } else {
            code_buf[code_buf_ptr++] = (unsigned char) match_position;
            code_buf[code_buf_ptr++] = (unsigned char)
                (((match_position >> 4) & 0xf0)
             | (match_length - (THRESHOLD + 1)));  /* Send position and
                    length pair. Note match_length > THRESHOLD. */
        }
        if ((mask <<= 1) == 0) {  /* Shift mask left one bit. */
            for (i = 0; i < code_buf_ptr; i++) {  /* Send at most 8 units of */
                _put(code_buf[i]);    /* code together */
            }
            codesize += code_buf_ptr;
            code_buf[0] = 0;  code_buf_ptr = mask = 1;
        }
        last_match_length = match_length;
        for (i = 0; i < last_match_length && ilen; i++) {
            _get(c);
            DeleteNode(lson, rson, dad, s);        /* Delete old strings and */
            text_buf[s] = c;    /* read new bytes */
            if (s < F - 1) text_buf[s + N] = c;  /* If the position is
                near the end of buffer, extend the buffer to make
                string comparison easier. */
            s = (s + 1) & (N - 1);  r = (r + 1) & (N - 1);
                /* Since this is a ring buffer, increment the position
                   modulo N. */
            InsertNode(text_buf, lson, rson, dad, r);    /* Register the string in text_buf[r..r+F-1] */
        }
        while (i++ < last_match_length) {    /* After the end of text, */
            DeleteNode(lson, rson, dad, s);                    /* no need to read, but */
            s = (s + 1) & (N - 1);  r = (r + 1) & (N - 1);
            if (--len) InsertNode(text_buf, lson, rson, dad, r);        /* buffer may not be empty. */
        }
    } while (len > 0);    /* until length of string to be processed is zero */
    if (code_buf_ptr > 1) {        /* Send remaining code. */
        for (i = 0; i < code_buf_ptr; i++) {
            _put(code_buf[i]);
        }
        codesize += code_buf_ptr;
    }

    return codesize;
}

#undef _put
#define _put(c) \
    *ostr++ = c;

size_t LzssDecode(const unsigned char* istr, size_t ilen, char *ostr, size_t olen)    /* Just the reverse of Encode(). */
{
    unsigned char text_buf[N + F - 1];          /* ring buffer of size N, with extra F-1 bytes to facilitate string comparison */
    int  i, j, k, r, c;
    unsigned int  flags;
    int limit = ilen;
    char *obak = ostr;
    
    for (i = 0; i < N - F; i++) text_buf[i] = 0;
    r = N - F;  flags = 0;
    for ( ; ; ) {
        if (((flags >>= 1) & 256) == 0) {
            _get(c);
            flags = c | 0xff00;        /* uses higher byte cleverly */
        }                            /* to count eight */
        if (flags & 1) {
            _get(c);
            _put(c);
            text_buf[r++] = c;  r &= (N - 1);
        } else {
            _get(i);
            _get(j);
            i |= ((j & 0xf0) << 4);  j = (j & 0x0f) + THRESHOLD;
            for (k = 0; k <= j; k++) {
                c = text_buf[(i + k) & (N - 1)];
                _put(c);
                text_buf[r++] = c;  r &= (N - 1);
            }
        }
    }
    return (size_t)(ostr - obak);
}

#undef _get
#undef _put

unsigned lzss_encode(const void *in, unsigned inlen, void *out, unsigned outlen, unsigned flags) {
    size_t rc = LzssEncode((const char*)in, (size_t)inlen, (char*)out, (size_t)outlen);
    return (unsigned)rc;
}
unsigned lzss_decode(const void *in, unsigned inlen, void *out, unsigned outlen) {
    size_t rc = LzssDecode((const unsigned char*)in, (size_t)inlen, (char*)out, (size_t)outlen);
    return (unsigned)rc;
}

#endif // LZSS_C


#ifdef LZSS_DEMO
#pragma once
#include <stdio.h>
int main() {
    const char *longcopy = "Hello world! Hello world! Hello world! Hello world!";

    int level=1;
    char out[128];
    size_t outlen = lzss_encode(longcopy, strlen(longcopy)+1, out, 128, level);
    printf("%s %d->%d\n", outlen ? "ok" : "fail", (int)strlen(longcopy)+1, (int)outlen);

    char redo[128];
    size_t unpacked = lzss_decode(out, outlen, redo, 128);
    printf("%d->%d %s\n", (int)outlen, (int)unpacked, redo);
}
#define main main__
#endif // LZSS_DEMO
