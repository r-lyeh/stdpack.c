# stdarc.c

This is a collection of public domain compressors from different authors. Sources have been minimally refurbished where due and/or backported from C++ to plain C. Most modifications were adding a standarized api convention and converting into single-header, converting from file to in-memory based, also added re-entrancy and thread-safety fixes, portability improvements, etc.

At this point, you do not likely want to use these compressors in your serious products. First-citizen compressors like Snappy, LZ4/LZ4HC, ZSTD, ZLIB, BZIP2 or LZMA SDK, are much more mature, more portable, extremely faster and way more battle-tested codebases than the ones found here. Anyways, if your main motivations are (like mine) having more liberal licensing terms, smaller codebases, ease of free standalone functions, or love header-only files or think they produce less hassles in general, etc. then you are welcome here.

By using these compressors you are not required to credit original authors neither to display licensing terms in your end-products, but I encourage you to acknowledge the original authors if possible, as they crafted with <3 these little gems over the years.

# About compression (tips)

1. When possible, feed each dedicated hardware its own compressed data format. Ie, feed BCn,ETCn,PVRTn... textures to the GPUs that support them, also feed AAC,MPn,OGG... bitstreams to the sound decoders that directly support them, and so on. No general compression algorithm can beat specific compression for specialized hardware. Always know your data and evolve your tools around it: as general rule, general compression should be always the last resort.

1. When picking compressors, ensure that memory consumption of your picked algorithm matches your target platform budget. For example, BALZ or LZMA with large dictionaries may take dozens or hundreds of MiBs while decompressing, so ensure they will run in a desktop and not in a gameboy in this case.

1. Once algorithms are chosen, evaluate both de/compression times. Think about main usage for each case. Servers and backend services could use the fastest compressor (to store data as fast as possible), game assets could use the fastest decompression (to minimize loading times), peer2peer data exchange could use both fastest de/compressors (to handle traffic quickly), and offline or installation data could definitely use the heaviest compressor (to minimize distribution costs), and so on.

1. Finally, consider transcoding to faster codecs if you plan to reload same data over and over. Ie, during your game installation you may transcode offline LZMA streams and save them as LZ4X for faster runtime loading.

1. TL;DR? Personal preferences? balz:1 for offline tasks, ulz:0 or lz4x:1 for runtime tasks, zlib:3 elsewhere.

# Some metrics (YMMV)

```
C:\prj\stdarc>demo.exe --0 \prj\enwik8
Y 100000000 ->    58769363 ppp  58.77% c:0.783s d:0.495s \prj\enwik8
Y 100000000 ->    49063365 lzss 49.06% c:14.17s d:0.424s \prj\enwik8
Y 100000000 ->    47873720 lzw3 47.87% c:1.143s d:0.372s \prj\enwik8
Y 100000000 ->    52660722 lz4x 52.66% c:1.004s d:0.165s \prj\enwik8
Y 100000000 ->    51940556 ulz  51.94% c:0.690s d:0.176s \prj\enwik8
Y 100000000 ->    63192823 zlib 63.19% c:1.837s d:0.541s \prj\enwik8
Y 100000000 ->    47639504 crsh 47.64% c:2.098s d:0.530s \prj\enwik8
Y 100000000 ->    30253237 balz 30.25% c:11.07s d:5.253s \prj\enwik8
Y 100000000 ->    39973598 lzma 39.97% c:13.65s d:2.810s \prj\enwik8
57.96s
```

```
C:\prj\stdarc>demo.exe --10 \prj\enwik8
Y 100000000 ->    58769363 ppp  58.77% c:0.783s d:0.474s \prj\enwik8
Y 100000000 ->    49063365 lzss 49.06% c:13.76s d:0.417s \prj\enwik8
Y 100000000 ->    47873720 lzw3 47.87% c:1.164s d:0.367s \prj\enwik8
Y 100000000 ->    45876420 lz4x 45.88% c:1.459s d:0.166s \prj\enwik8
Y 100000000 ->    40285775 ulz  40.29% c:11.95s d:0.169s \prj\enwik8
Y 100000000 ->    36465845 zlib 36.47% c:6.581s d:0.462s \prj\enwik8
Y 100000000 ->    31856827 crsh 31.86% c:333.2s d:0.419s \prj\enwik8
Y 100000000 ->    28460117 balz 28.46% c:32.36s d:4.989s \prj\enwik8
Y 100000000 ->    26702913 lzma 26.70% c:86.37s d:1.577s \prj\enwik8
497.37s
```

# Project goals
- C.
- Drop & use.
- Small codebase.
- Good enough performance.
- Public domain or unlicensed source code. No attribution required.

# Todo
- Test other environments. Currently VS2019 only.
- Add small compression filters (like E8E9).
- Add small file archivers (pak, zip, tar and custom). Maybe.
- Add vfs and mounting support. Maybe.
- Optimize LZSS. Could compressor be nearly as performant as LZJB with some effort?
- Optimize BALZ. Could those in-mem bufferings be skipped?
- Consider ZPAQ (large C++), BCM (divsufsort is MIT licensed) and yalz77 too. LZP1, LZJB maybe?
- Bugfixes and cleaning up.
