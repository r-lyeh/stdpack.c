// directory iteration.
// - rlyeh, public domain.

#ifndef DIR_H
#define DIR_H

typedef struct dir dir;

dir *dir_open(const char *filename, const char *mode); // recursive 'r'

    int dir_find(dir*, const char *entryname); // returns entry number; or <0 if not found.
    unsigned dir_count(dir*);
        char*    dir_name(dir*, unsigned index);
        unsigned dir_size(dir*, unsigned index);
        unsigned dir_file(dir*, unsigned index); // dir_isfile? bool?
        void*    dir_read(dir*, unsigned index); // must free() after use

void dir_close(dir*);

#endif

// -----------------------------------------------------------------------------

#ifdef DIR_C
#pragma once
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#  if defined _WIN32 && defined(__TINYC__)
#include <windows.h>  // tcc
#elif defined _WIN32
#include <winsock2.h> // msc+gcc
#else
#include <dirent.h>
#endif

#ifndef STRDUP
#define STRDUP strdup
#endif

#ifndef REALLOC
#define REALLOC realloc
#endif

#ifndef ERR
#define ERR(NUM, ...) (fprintf(stderr, "" __VA_ARGS__), fprintf(stderr, "(%s:%d)\n", __FILE__, __LINE__), fflush(stderr), (NUM)) // (NUM)
#endif

typedef struct dir_entry {
    char *filename;
    size_t size;
    size_t is_dir : 1;
} dir_entry;

struct dir {
    dir_entry *entry;
    unsigned count;
};

// ---

#if !defined(S_ISDIR)
#   define S_ISDIR(mode) (((mode) & S_IFMT) == S_IFDIR)
#endif

int dir_yield(dir *d, const char *pathfile, char *name, int namelen) {
    int ok = 0;
#ifdef _WIN32
    WIN32_FIND_DATAA fdata = { 0 };
    snprintf(name, namelen, "%s/*", pathfile);
    for( HANDLE h = FindFirstFileA(name, &fdata ); h != INVALID_HANDLE_VALUE; (ok = FindClose( h ), h = INVALID_HANDLE_VALUE, 1)) {
        for( int next = 1; next; next = FindNextFileA(h, &fdata) != 0 ) {
            if( fdata.cFileName[0] == '.' ) continue;
            int is_dir = (fdata.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) > 0;
            snprintf(name, namelen, "%s/%s%s", pathfile, fdata.cFileName, is_dir ? "/" : "");
            struct stat st; if( !is_dir ) if(stat(name, &st) < 0) continue;
            // add
            dir_entry de = { STRDUP(name), is_dir ? 0 : st.st_size, is_dir };
            d->entry = (dir_entry*)REALLOC(d->entry, ++d->count * sizeof(dir_entry));
            d->entry[d->count-1] = de;
        }
    }
#else
    snprintf(name, namelen, "%s/", pathfile);
    for( DIR *dir = opendir(name); dir; ok = (closedir(dir), dir = 0, 1)) {
        for( struct dirent *ep; ep = readdir(dir); ) {
            if( ep->d_name[0] == '.' ) continue;
            snprintf(name, namelen, "%s/%s", pathfile, ep->d_name);
            struct stat st; if( stat(name, &st) < 0 ) continue;
            DIR *tmp = opendir(ep->d_name); int is_dir = !!tmp; if(tmp) closedir(tmp);
            // if( is_dir && recursive ) { dir_yield(d,name); }
            // add
            dir_entry de = { STRDUP(name), is_dir ? 0 : st.st_size, is_dir };
            d->entry = (dir_entry*)REALLOC(d->entry, ++d->count * sizeof(dir_entry));
            d->entry[d->count-1] = de;
        }
    }
#endif
    return ok;
}

dir *dir_open(const char *pathfile, const char *mode) {
    // if(mode[0] == 'R') return ERR(NULL, "(R) not supported for now");
    dir *d = (dir*)REALLOC(0, sizeof(dir)), zero = {0}; *d = zero;

    char *clean = STRDUP( pathfile );
    for( int i = 0; clean[i]; ++i ) if(clean[i] == '\\') clean[i] = '/';
    for( int len = strlen(clean); clean[--len] == '/'; ) clean[len] = '\0';

    char buffer[2048];
    dir_yield(d, clean, buffer, 2048);

    REALLOC(clean, 0);
    return d;
}

int dir_find(dir *d, const char *entryname) {
    for( int i = d->count; --i >= 0; ) { // in case of several copies, grab most recent file (last coincidence)
        if( 0 == strcmp(entryname, d->entry[i].filename)) return i;
    }
    return -1;
}

unsigned dir_count(dir *d) {
    return d ? d->count : 0;
}

char* dir_name(dir *d, unsigned index) {
    return d && index < d->count ? d->entry[index].filename : 0;
}

unsigned dir_size(dir *d, unsigned index) {
    return d && index < d->count ? (unsigned)d->entry[index].size : 0;
}

unsigned dir_file(dir *d, unsigned index) {
    return d && index < d->count ? (unsigned)!d->entry[index].is_dir : 0;
}

void *dir_read(dir *d, unsigned index) {
    if( d && index < d->count ) {
        void *data = 0;
        for( FILE *fp = fopen(d->entry[index].filename, "rb"); fp; fclose(fp), fp = 0) {
            size_t len = d->entry[index].size;
            data = REALLOC(0, len);
            if( data && fread(data, 1, len, fp) != len ) {
                data = REALLOC(data, 0);
            }
        }
        return data;
    }
    return 0;
}

void dir_close(dir *d) {
    for( int i = 0; i < d->count; ++i) {
        REALLOC(d->entry[i].filename, 0);
    }
    dir zero = {0};
    *d = zero;
    REALLOC(d, 0);
}

#ifdef DIR_DEMO
int main( int argc, char **argv ) {
    dir *d = dir_open(argc > 1 ? argv[1] : "./", "rb");
    if( d ) {
        for( int i = 0; i < dir_count(d); ++i ) {
            if( dir_file(d,i) )
            printf("%3d) %11d %s\n", i + 1, dir_size(d,i), dir_name(d,i));
            else
            printf("%3d) %11s %s\n", i + 1, "<dir>", dir_name(d,i));
            char *data = dir_read(d,i);
            if(argc > 2 && !strcmp(argv[2],dir_name(d,i))) printf("%.*s\n", dir_size(d,i), data);
            free(data);
        }
        dir_close(d);
    }
}
#define main main__
#endif //DIR_DEMO
#endif //DIR_C
