// gnu tar and ustar extraction
// - rlyeh, public domain.

#ifndef TAR_H
#define TAR_H

typedef struct tar tar;

tar *tar_open(const char *filename, const char *mode);

    int tar_find(tar*, const char *entryname); // returns entry number; or <0 if not found.
    unsigned tar_count(tar*);
        char*    tar_name(tar*, unsigned index);
        unsigned tar_size(tar*, unsigned index);
        unsigned tar_offset(tar*, unsigned index);
        void*    tar_extract(tar*, unsigned index); // must free() after use

void tar_close(tar *t);

#endif

// -----------------------------------------------------------------------------

#ifdef TAR_C
#pragma once
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifndef STRDUP
#define STRDUP strdup
#endif

#ifndef REALLOC
#define REALLOC realloc
#endif

#ifndef ERR
#define ERR(NUM, ...) (fprintf(stderr, "" __VA_ARGS__), fprintf(stderr, "(%s:%d)\n", __FILE__, __LINE__), fflush(stderr), (NUM)) // (NUM)
#endif

struct tar {
    FILE *in;
    unsigned count;
    struct tar_entry {
    char *filename;
    unsigned size;
    size_t offset;
    } *entries;
};

// equivalent to sscanf(buf, 8, "%.7o", &size); or (12, "%.11o", &modtime)
// ignores everything after first null or space, including trailing bytes
uint64_t tar__octal( const char *src, const char *eof ) {
    uint64_t sum = 0, mul = 1;
    const char *ptr = eof;
    while( ptr-- >= src ) eof  = ( 0 != ptr[1] && 32 != ptr[1] ) ? eof : ptr;
    while( eof-- >= src ) sum += (uint8_t)(eof[1] - '0') * mul, mul *= 8;
    return sum;
}

typedef int (*tar_callback)(const char *filename, unsigned inlen, size_t offset, void *userdata);

int tar__push_entry(const char *filename, unsigned inlen, size_t offset, void *userdata) {
    tar *t = (tar *)userdata;

    unsigned index = t->count;
    t->entries = REALLOC(t->entries, (++t->count) * sizeof(struct tar_entry));
    struct tar_entry *e = &t->entries[index];

    e->filename = STRDUP(filename);
    e->size = inlen;
    e->offset = offset;

    return 1;
}

int tar__parse( FILE *in, tar_callback yield, void *userdata ) {
    enum {
        name     =   0, // (null terminated)
        mode     = 100, // (octal)
        uid      = 108, // (octal)
        gid      = 116, // (octal)
        size     = 124, // (octal)
        modtime  = 136, // (octal)
        checksum = 148, // (octal)
        type     = 156, // \0|'0':file,1:hardlink,2:symlink,3:chardev,4:blockdev,5:dir,6:fifo,L:longnameblocknext
        linkname = 157, // if !ustar link indicator
        magic    = 257, // if ustar "ustar" -- 6th character may be space or null, else zero
        version  = 263, // if ustar "00", else zero
        uname    = 265, // if ustar owner username, else zero
        gname    = 297, // if ustar owner groupname, else zero
        devmajor = 329, // if ustar device major number, else zero
        devminor = 337, // if ustar device minor number , else zero
        path     = 345, // if ustar filename prefix, else zero
        padding  = 500, // if ustar relevant for checksum, else zero
        total    = 512
    };
    // handle both regular tar and ustar tar filenames until end of tar is found
    char header[512], entry[512], blank[512] = {0};
    while( !ferror(in) ) {
        if( 512 != fread(header, 1, 512, in ) ) break;
        if( memcmp( header, blank, 512 ) ) {                                      // if not end of tar
            if( !memcmp( header+magic, "ustar", 5 ) ) {                           // if valid ustar
                int namelen = strlen(header+name), pathlen = strlen(header+path); // read filename
                snprintf(entry, 512, "%.*s" "%s" "%.*s",
                    pathlen < 155 ? pathlen : 155, header+path,
                    pathlen ? "/" : "",
                    namelen < 100 ? namelen : 100, header+name );
                switch( header[type] ) {
                    default:                                                      // unsupported file type
                    break; case '5': //yield(entry.back()!='/'?entry+'/':entry,0);// directory
                    break; case 'L': strcpy(entry, header+name); fread(header,1,512,in); // gnu tar long filename
                    break; case '0': case 0: {                                    // regular file
                        uint64_t len = tar__octal(header+size, header+modtime);    // decode octal size
                        int cont = yield(entry, len, ftell(in), userdata);        // yield entry
                        fseek(in,len,SEEK_CUR);                                   // skip blob
                        fseek(in,(512 - (len & 511)) & 511,SEEK_CUR);             // skip padding
                    }
                }
            } else return ERR(0, "not a .tar file");
        } else return ferror(in) ? ERR(0, "file error") : 1;
    }
    return ERR(0, "read error");
}

// ---

tar *tar_open(const char *filename, const char *mode) {
    if(mode[0] != 'r') return ERR(NULL, "(w) and (a) not supported for now");
    FILE *in = fopen(filename, "rb");
    if(!in) return ERR(NULL, "cant open file '%s' for reading", filename);

    tar zero = {0}, *t = REALLOC(0, sizeof(tar));
    if( !t ) { fclose(in); return ERR(NULL, "out of mem"); }

    *t = zero;
    t->in = in;
    tar__parse(in, tar__push_entry, t);
    return t;
}

int tar_find(tar *t, const char *entryname) {
    if( t->in ) for( int i = t->count; --i >= 0; ) { // in case of several copies, grab most recent file (last coincidence)
        if( 0 == strcmp(entryname, t->entries[i].filename)) return i;
    }
    return -1;
}

unsigned tar_count(tar *t) {
    return t ? t->count : 0;
}

char* tar_name(tar *t, unsigned index) {
    return t && index < t->count ? t->entries[index].filename : 0;
}

unsigned tar_size(tar *t, unsigned index) {
    return t && index < t->count ? t->entries[index].size : 0;
}

unsigned tar_offset(tar *t, unsigned index) {
    return t && index < t->count ? (unsigned)t->entries[index].offset : 0;
}

void *tar_extract(tar *t, unsigned index) {
    if( t && index < t->count ) {
        fseek(t->in, t->entries[index].offset, SEEK_SET);
        size_t len = t->entries[index].size;
        void *data = REALLOC(0, t->entries[index].size);
        fread(data, 1, len, t->in);
        return data;
    }
    return 0;
}

void tar_close(tar *t) {
    fclose(t->in);
    for( int i = 0; i < t->count; ++i) {
        REALLOC(t->entries[i].filename, 0);
    }
    tar zero = {0};
    *t = zero;
    REALLOC(t, 0);
}

#ifdef TAR_DEMO
int main( int argc, char **argv ) {
    if(argc <= 1) exit(printf("%s file.tar [file_to_view]\n", argv[0]));
    tar *t = tar_open(argv[1], "rb");
    if( t ) {
        for( int i = 0; i < tar_count(t); ++i ) {
            printf("%d) %s (%u bytes)\n", i+1, tar_name(t,i), tar_size(t,i));
            char *data = tar_extract(t,i);
            if(argc>2) if(0==strcmp(argv[2],tar_name(t,i))) printf("%.*s\n", tar_size(t,i), data);
            free(data);
        }
        tar_close(t);
    }
}
#define main main__
#endif //TAR_DEMO
#endif //TAR_C
