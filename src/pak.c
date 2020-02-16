// pak file reading/writing/appending.
// - rlyeh, public domain.
//
// ## PAK
// - [ref] https://quakewiki.org/wiki/.pak (public domain).
// - Header: 12 bytes
//   - "PACK"           4-byte
//   - directory offset uint32
//   - directory size   uint32 (number of files by dividing this by 64, which is sizeof(pak_entry))
//
// - File Directory Entry (Num files * 64 bytes)
//   - Each Directory Entry: 64 bytes
//     - file name     56-byte null-terminated string. Includes path. Example: "maps/e1m1.bsp".
//     - file offset   uint32 from beginning of pak file.
//     - file size     uint32

#ifndef PAK_H
#define PAK_H

typedef struct pak pak;

pak* pak_open(const char *fname, const char *mode /*a,r,w*/);

    // (w)rite or (a)ppend modes only
    int pak_append_file(pak*, const char *filename, FILE *in);
    int pak_append_data(pak*, const char *filename, const void *in, unsigned inlen);
    
    // (r)ead only mode
    int pak_find(pak*,const char *fname); // return <0 if error; index otherwise.
    unsigned pak_count(pak*);
        unsigned pak_size(pak*,unsigned index);
        unsigned pak_offset(pak*, unsigned index);
        char *pak_name(pak*,unsigned index);
        void *pak_extract(pak*, unsigned index); // must free() after use

void pak_close(pak*);

#endif

// ---

#ifdef PAK_C
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifndef REALLOC
#define REALLOC realloc
#endif

#ifndef ERR
#define ERR(NUM, ...) (fprintf(stderr, "" __VA_ARGS__), fprintf(stderr, "(%s:%d)\n", __FILE__, __LINE__), fflush(stderr), (NUM)) // (NUM)
#endif

#include <stdint.h>
static inline uint32_t pak_swap32( uint32_t t ) { return (t >> 24) | (t << 24) | ((t >> 8) & 0xff00) | ((t & 0xff00) << 8); }

#if defined(_M_IX86) || defined(_M_X64) // #ifdef LITTLE
#define htob32(x) pak_swap32(x)
#define btoh32(x) pak_swap32(x)
#define htol32(x) (x)
#define ltoh32(x) (x)
#else
#define htob32(x) (x)
#define btoh32(x) (x)
#define htol32(x) pak_swap32(x)
#define ltoh32(x) pak_swap32(x)
#endif

#pragma pack(push, 1)

typedef struct pak_header {
    char id[4];
    uint32_t offset;
    uint32_t size;
} pak_header;

typedef struct pak_file {
    char name[56];
    uint32_t offset;
    uint32_t size;
} pak_file;

#pragma pack(pop)

typedef int static_assert_sizeof_pak_header[sizeof(pak_header) == 12];
typedef int static_assert_sizeof_pak_file[sizeof(pak_file) == 64];

typedef struct pak {
    FILE *in, *out;
    int dummy;
    pak_file *entries;
    unsigned count;
} pak;

pak *pak_open(const char *fname, const char *mode) {
    struct stat buffer;
    int exists = (stat(fname, &buffer) == 0);
    if(mode[0] == 'a' && !exists ) mode = "wb";

    if(mode[0] != 'w' && mode[0] != 'r' && mode[0] != 'a') return NULL;

    FILE *fp = fopen(fname, mode[0] == 'w' ? "wb" : mode[0] == 'r' ? "rb" : "r+b");
    if(!fp) return ERR(NULL, "cant open file '%s' in '%s' mode", fname, mode);

    pak *p = malloc(sizeof(pak)), zero = {0};
    if(!p) return fclose(fp), ERR(NULL, "out of mem");
    *p = zero;

    if( mode[0] == 'r' || mode[0] == 'a' ) {
        pak_header header = {0};

        if( fread(&header, 1, sizeof(pak_header), fp) != sizeof(pak_header) ) {
            return fclose(fp), ERR(NULL, "read error");
        }
        if( memcmp(header.id, "PACK", 4) ) {
            return fclose(fp), ERR(NULL, "not a .pak file");
        }

        header.offset = ltoh32(header.offset);
        header.size = ltoh32(header.size);

        unsigned num_files = header.size / sizeof(pak_file);

        if( fseek(fp, header.offset, SEEK_SET) != 0 ) {
            return fclose(fp), ERR(NULL, "read error");
        }

        p->count = num_files;
        p->entries = REALLOC(0, num_files * sizeof(pak_file));

        if( fread(p->entries, num_files, sizeof(pak_file), fp) != sizeof(pak_file) ) {
            goto fail;
        }

        for( unsigned i = 0; i < num_files; ++i ) {
            pak_file *e = &p->entries[i];
            e->offset = ltoh32(e->offset);
            e->size = ltoh32(e->size);
        }

        if( mode[0] == 'a' ) {
            // resize (by truncation)
            size_t resize = header.offset;
            int fd = fileno(fp);
            if( fd != -1 ) {
                #ifdef _WIN32
                    int ok = 0 == _chsize_s( fd, resize );
                #else
                    int ok = 0 == ftruncate( fd, (off_t)resize );
                #endif
                fflush(fp);
                fseek( fp, 0L, SEEK_END );
            }

            p->out = fp;
            p->in = NULL;
        } else {
            p->in = fp;
        }

        return p;
    }


    if(mode[0] == 'w') {
        p->out = fp;

        // write temporary header
        char header[12] = {0};
        if( fwrite(header, 1,12, p->out) != 12) goto fail;

        return p;
    }

fail:;
    if(fp) fclose(fp);
    if(p->entries) REALLOC(p->entries, 0);
    if(p) REALLOC(p, 0);

    return NULL;
}

int pak_append_data(pak *p, const char *filename, const void *in, unsigned inlen) {
    if(!p->out) return ERR(0, "read-only pak file");

    // index meta
    unsigned index = p->count++;
    p->entries = REALLOC(p->entries, p->count * sizeof(pak_file));
    pak_file *e = &p->entries[index], zero = {0};
    *e = zero;
    snprintf(e->name, 55, "%s", filename); // @todo: verify 56 chars limit
    e->size = inlen;
    e->offset = ftell(p->out);

    // write blob
    fwrite(in, 1, inlen, p->out);

    return !ferror(p->out);
}

int pak_append_file(pak *p, const char *filename, FILE *in) {
    // index meta
    unsigned index = p->count++;
    p->entries = REALLOC(p->entries, p->count * sizeof(pak_file));
    pak_file *e = &p->entries[index], zero = {0};
    *e = zero;
    snprintf(e->name, 55, "%s", filename); // @todo: verify 56 chars limit
    e->offset = ftell(p->out);

    char buf[1<<15];
    while(!feof(in) && !ferror(in)) {
        size_t bytes = fread(buf, 1, sizeof(buf), in);
        fwrite(buf, 1, bytes, p->out);
    }

    e->size = ftell(p->out) - e->offset;

    return !ferror(p->out);
}


void pak_close(pak *p) {
    if(p->out) {
        // write toc
        uint32_t seek = 0 + 12, dirpos = (uint32_t)ftell(p->out), dirlen = p->count * 64;
        for(unsigned i = 0; i < p->count; ++i) {
            pak_file *e = &p->entries[i];
            // write name (truncated if needed), and trailing zeros
            char zero[56] = {0};
            int namelen = strlen(e->name);
            fwrite( e->name, 1, namelen >= 56 ? 55 : namelen, p->out );
            fwrite( zero, 1, namelen >= 56 ? 1 : 56 - namelen, p->out );
            // write offset + length pair
            uint32_t pseek = htol32(seek);    fwrite( &pseek, 1,4, p->out );
            uint32_t psize = htol32(e->size); fwrite( &psize, 1,4, p->out );
            seek += e->size;
        }

        // patch header
        fseek(p->out, 0L, SEEK_SET);
        fwrite("PACK", 1,4, p->out);
        dirpos = htol32(dirpos); fwrite( &dirpos, 1,4, p->out );
        dirlen = htol32(dirlen); fwrite( &dirlen, 1,4, p->out );
    }

    // close streams
    if(p->in) fclose(p->in);
    if(p->out) fclose(p->out);

    // clean up
    for(unsigned i = 0; i < p->count; ++i) {
        pak_file *e = &p->entries[i];
    }
    REALLOC(p->entries, 0);

    // delete
    pak zero = {0};
    *p = zero;
    REALLOC(p, 0);
}

int pak_find(pak *p, const char *filename) {
    if( p->in ) {
        for( int i = p->count; --i >= 0; ) {
            if(!strcmp(p->entries[i].name, filename)) return i;
        }
    }
    return -1;
}

unsigned pak_count(pak *p) {
    return p->in ? p->count : 0;
}

unsigned pak_offset(pak *p, unsigned index) {
    return p->in && index < p->count ? p->entries[index].offset : 0;
}

unsigned pak_size(pak *p, unsigned index) {
    return p->in && index < p->count ? p->entries[index].size : 0;
}

char *pak_name(pak *p, unsigned index) {
    return p->in && index < p->count ? p->entries[index].name : NULL;
}

void *pak_extract(pak *p, unsigned index) {
    if( p->in && index < p->count ) {
        pak_file *e = &p->entries[index];
        if( fseek(p->in, e->offset, SEEK_SET) != 0 ) {
            return ERR(NULL, "cant seek");
        }
        void *buffer = REALLOC(0, e->size);
        if( !buffer ) {
            return ERR(NULL, "out of mem");
        }
        if( fread(buffer, 1, e->size, p->in) != e->size ) {
            REALLOC(buffer, 0);
            return ERR(NULL, "cant read");
        }
        return buffer;
    }
    return NULL;
}

#ifdef PAK_DEMO
int main(int argc, char **argv) {
    puts("creating test.pak archive (3) ...");
    pak *p = pak_open("test.pak", "wb");
    if( p ) {
        pak_append_data(p, "/index.txt", "just a file", strlen("just a file"));
        pak_append_data(p, "/file/name1.txt", "just another file #1", strlen("just another file #1"));
        pak_append_data(p, "/file/name2.txt", "just another file #2", strlen("just another file #2"));
        pak_close(p);
    }

    puts("appending file to test.pak (1) ...");
    p = pak_open("test.pak", "a+b");
    if( p ) {
        pak_append_data(p, "/new/file", "this is an appended file", strlen("this is an appended file"));
        pak_close(p);
    }

    const char *fname = argc > 1 ? argv[1] : "test.pak";
    printf("listing %s archive ...\n", fname);
    p = pak_open(fname, "rb");
    if( p ) {
        for( unsigned i = 0; i < pak_count(p); ++i ) {
            printf("  %d) @%08x %11u %s ", i+1, pak_offset(p,i), pak_size(p,i), pak_name(p,i));
            void *data = pak_extract(p,i);
            printf("\r%c\n", data ? 'Y':'N');
            if(argc > 2 && data) 
                if(i == pak_find(p,argv[2]))
                    printf("%.*s\n", (int)pak_size(p,i), (char*)data);
            free(data);
        }
        pak_close(p);
    }

    puts("ok");
}
#endif // PAK_DEMO
#endif // PAK_C
