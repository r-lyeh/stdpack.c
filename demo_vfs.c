#include <stdio.h>

#define STDARC_C
#include "stdarc.c"

int main() {
    vfs_mount("src/"); // directories/must/end/with/slash/
    vfs_mount("demo.zip"); // zips supported
    printf("vfs.c file found? %s\n", vfs_load("vfs.c", 0) ? "Y":"N");
    printf("stdarc.c file found? %s\n", vfs_load("stdarc.c", 0) ? "Y":"N");
    printf("demo_zip.c file found? %s\n", vfs_load("demo_zip.c", 0) ? "Y":"N");
}
