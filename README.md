# stdpack.c

This is [a collection of public domain compressors](src) from different authors. Sources have been minimally refurbished where due and/or backported from C++ to plain C. Most modifications were adding a standarized api convention and converting into single-header, converting from file to in-memory based, also added re-entrancy and thread-safety fixes, portability improvements, etc.

At this point, you do not likely want to use these compressors in serious projects. First-citizen compressors like Snappy, LZ4/LZ4HC, ZSTD, ZLIB, BZIP2 or LZMA SDK are generally better bets: way more faster, mature, portable and battle-tested than the snippets found here. Still, you are welcome to use these snippets if your motivations are similar to mine: to have less headaches by embedding 3rd party software with no clauses, by preferring header-only files, by not hand-crafting makefiles, by avoid dealing with dynamic libraries, and by inspecting and debugging small codebases in general.

By using these compressors you are not required to credit original authors neither to display licensing terms in your final products, but I encourage you to acknowledge the original authors if possible, as they crafted with <3 these little gems over the years.

# Archives (zip, tar, pak, ...)

Please check https://github.com/r-lyeh/stdarc.c instead

# About compression (tips)

1. When possible, feed each dedicated hardware its own compressed data format. Ie, feed BCn,ETCn,PVRTn... textures to the GPUs that support them, also feed AAC,MPn,OGG... bitstreams to the sound decoders that directly support them, and so on. No general compression algorithm can beat specific compression for specialized hardware. Always know your data and evolve your tools around it. As general rule, general compression should be always the last resort.

1. When picking compressors, know your algorithms and your platform budgets. Ensure that memory consumption of your picked algorithm matches your target platform budget. For example, BALZ or LZMA with large dictionaries will take dozens or hundreds of MiBs while decompressing, so ensure they will be used in a desktop computer and not in a 512M mobile.

1. Once algorithms are chosen, evaluate both de/compression times. Think about main usage for each case. Servers and backend services could use the fastest compressor (to store data as fast as possible), game assets could use the fastest decompression (to minimize loading times), peer2peer data exchange could use both fastest de/compressors (to handle traffic quickly), and offline or installation data could definitely use the heaviest compressor (to minimize distribution costs), and so on.

1. Finally, consider transcoding to faster codecs if you plan to reload same data over and over. Ie, during your game installation you may transcode offline LZMA streams and save them as LZ4X for faster runtime loading.

1. TL;DR? Personal preferences? balz:1 for offline tasks, ulz:0 or lz4x:1 for runtime tasks, defl:3 elsewhere.

# Some metrics (YMMV)

```
C:\prj\stdpack>test.exe --benchmark --0 \prj\enwik8 | sort
Y 100000000 ->    58555380  ppp.0 58.56% c:0.721s d:0.413s \prj\enwik8
Y 100000000 ->    56013654 lzp1.0 56.01% c:0.821s d:0.537s \prj\enwik8
Y 100000000 ->    52650827 lz4x.0 52.65% c:0.754s d:0.097s \prj\enwik8
Y 100000000 ->    51923215  ulz.0 51.92% c:0.646s d:0.094s \prj\enwik8
Y 100000000 ->    49061221 lzss.0 49.06% c:14.24s d:0.355s \prj\enwik8
Y 100000000 ->    47874204 lzw3.0 47.87% c:1.144s d:0.315s \prj\enwik8
Y 100000000 ->    47620393 crsh.0 47.62% c:1.948s d:0.468s \prj\enwik8
Y 100000000 ->    39963844 lzma.0 39.96% c:14.14s d:2.788s \prj\enwik8
Y 100000000 ->    39879472 defl.0 39.88% c:1.567s d:0.623s \prj\enwik8
Y 100000000 ->    30056093 balz.0 30.06% c:10.80s d:5.403s \prj\enwik8
Y 100000000 ->    27700930  bcm.0 27.70% c:14.15s d:12.36s \prj\enwik8
84.35s
```

```
C:\prj\stdpack>test.exe --benchmark --9 \prj\enwik8 | sort
Y 100000000 ->    58555380  ppp.9 58.56% c:0.733s d:0.419s \prj\enwik8
Y 100000000 ->    56013654 lzp1.9 56.01% c:0.819s d:0.532s \prj\enwik8
Y 100000000 ->    49061221 lzss.9 49.06% c:14.39s d:0.349s \prj\enwik8
Y 100000000 ->    47874204 lzw3.9 47.87% c:1.112s d:0.306s \prj\enwik8
Y 100000000 ->    46021464 lz4x.9 46.02% c:1.400s d:0.105s \prj\enwik8
Y 100000000 ->    40253295  ulz.9 40.25% c:12.71s d:0.104s \prj\enwik8
Y 100000000 ->    36407448 defl.9 36.41% c:4.373s d:0.559s \prj\enwik8
Y 100000000 ->    32038271 crsh.9 32.04% c:116.9s d:0.381s \prj\enwik8
Y 100000000 ->    28232820 balz.9 28.23% c:34.84s d:5.221s \prj\enwik8
Y 100000000 ->    25801721 lzma.9 25.80% c:106.6s d:1.516s \prj\enwik8
Y 100000000 ->    20789660  bcm.9 20.79% c:20.79s d:19.95s \prj\enwik8
344.09s
```

# Features
- De/compress memory blobs, with different compression algorithms.
- De/compress files, with different compressors algorithms.
- De/compress files in chunks, with different compression algorithms per chunk.

# Project goals
- C.
- Drop & use.
- Small codebase.
- Good enough performance.
- Public domain or unlicensed source code. No attribution required.

# Todo
- Test other environments. Currently VS2019+GCC+Clang only.
- Add small compression filters (like E8E9).
- Optimize LZSS. Could compressor be nearly as performant as LZJB with some effort?
- Optimize BALZ. Could those in-mem bufferings be skipped?
- Consider ZPAQ (large C++), ~~BCM~~, ~~yalz77~~. ~~LZJB maybe?~~
- Bugfixes and cleaning up.

# Links
- Ilya Muravyov: bcm, balz, crush, ulz, lz4x - https://github.com/encode84
- Micha Mettke: sinfl, sdefl - https://github.com/vurtun/lib
- Igor Pavlov: lzma - https://www.7-zip.org/
- Charles Bloom: lzp1 - https://www.cbloom.com/src/index_lz.html
- Ross Williams: lzrw3a - http://ross.net/compression/lzrw3a.html
- Haruhiko Okumura: lzss - https://oku.edu.mie-u.ac.jp/~okumura/compression/
- Dave Rand: ppp - https://tools.ietf.org/html/rfc1978
