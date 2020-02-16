// virtual filesystem (registered directories and/or compressed zip archives).
// - rlyeh, public domain.
//
// - note: vfs_mount() order matters (the most recent the higher priority).

void  vfs_mount(const char *path); // zipfile or directory/with/trailing/slash/
char* vfs_load(const char *filename, int *size); // must free() after use

// -----------------------------------------------------------------------------

#ifdef VFS_C
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static
char *vfs_read_file(const char *filename, int *len) {
    FILE *fp = fopen(filename, "rb");
    if( fp ) {
        fseek(fp, 0L, SEEK_END);
        size_t sz = ftell(fp);
        fseek(fp, 0L, SEEK_SET);
        char *bin = REALLOC(0, sz+1);
        fread(bin,sz,1,fp);
        fclose(fp);
        bin[sz] = 0;
        if(len) *len = (int)sz;
        return bin;
    }
    return 0;
}

typedef struct vfs_dir {
    char* path;
    // const 
    zip* archive;
    int is_directory;
    struct vfs_dir *next;
} vfs_dir;

static vfs_dir *dir_head = NULL;

void vfs_mount(const char *path) {
    zip *z = NULL;
    int is_directory = ('/' == path[strlen(path)-1]);
    if( !is_directory ) z = zip_open(path, "rb");
    if( !is_directory && !z ) return;

    vfs_dir *prev = dir_head, zero = {0};
    *(dir_head = REALLOC(0, sizeof(vfs_dir))) = zero;
    dir_head->next = prev;
    dir_head->path = STRDUP(path);
    dir_head->archive = z;
    dir_head->is_directory = is_directory;
}

char *vfs_load(const char *filename, int *size) { // must free() after use
    char *data = NULL;
    for(vfs_dir *dir = dir_head; dir && !data; dir = dir->next) {
        if( dir->is_directory ) {
            char buf[512];
            snprintf(buf, sizeof(buf), "%s%s", dir->path, filename);
            data = vfs_read_file(buf, size);
        } else {
            int index = zip_find(dir->archive, filename);
            data = zip_extract(dir->archive, index);
            if( size ) *size = zip_size(dir->archive, index);
        }
        // printf("%c trying %s in %s ...\n", data ? 'Y':'N', filename, dir->path);
    }
    return data;
}

#ifdef VFS_DEMO
int main() {
    vfs_mount("src/"); // directories/must/end/with/slash/
    vfs_mount("demo.zip"); // zips supported
    printf("vfs.c file found? %s\n", vfs_load("vfs.c", 0) ? "Y":"N");
    printf("stdarc.c file found? %s\n", vfs_load("stdarc.c", 0) ? "Y":"N");
    printf("demo_zip.c file found? %s\n", vfs_load("demo_zip.c", 0) ? "Y":"N");
}
#define main main__
#endif // VFS_DEMO
#endif // VFS_C
