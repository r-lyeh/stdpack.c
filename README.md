# stdarc.c

This is [a collection of public domain compressors](src) from different authors. Sources have been minimally refurbished where due and/or backported from C++ to plain C. Most modifications were adding a standarized api convention and converting into single-header, converting from file to in-memory based, also added re-entrancy and thread-safety fixes, portability improvements, etc.

At this point, you do not likely want to use these compressors in serious projects. First-citizen compressors like Snappy, LZ4/LZ4HC, ZSTD, ZLIB, BZIP2 or LZMA SDK are generally better bets: way more faster, mature, portable and battle-tested than the snippets found here. Still, you are welcome to use these snippets if your motivations are similar to mine: to have less headaches by embedding 3rd party software with no clauses, by preferring header-only files, by not hand-crafting makefiles, by avoid dealing with dynamic libraries, and by inspecting and debugging small codebases in general.

By using these compressors you are not required to credit original authors neither to display licensing terms in your final products, but I encourage you to acknowledge the original authors if possible, as they crafted with <3 these little gems over the years.

# About compression (tips)

1. When possible, feed each dedicated hardware its own compressed data format. Ie, feed BCn,ETCn,PVRTn... textures to the GPUs that support them, also feed AAC,MPn,OGG... bitstreams to the sound decoders that directly support them, and so on. No general compression algorithm can beat specific compression for specialized hardware. Always know your data and evolve your tools around it. As general rule, general compression should be always the last resort.

1. When picking compressors, know your algorithms and your platform budgets. Ensure that memory consumption of your picked algorithm matches your target platform budget. For example, BALZ or LZMA with large dictionaries will take dozens or hundreds of MiBs while decompressing, so ensure they will be used in a desktop computer and not in a 512M mobile.

1. Once algorithms are chosen, evaluate both de/compression times. Think about main usage for each case. Servers and backend services could use the fastest compressor (to store data as fast as possible), game assets could use the fastest decompression (to minimize loading times), peer2peer data exchange could use both fastest de/compressors (to handle traffic quickly), and offline or installation data could definitely use the heaviest compressor (to minimize distribution costs), and so on.

1. Finally, consider transcoding to faster codecs if you plan to reload same data over and over. Ie, during your game installation you may transcode offline LZMA streams and save them as LZ4X for faster runtime loading.

1. TL;DR? Personal preferences? balz:1 for offline tasks, ulz:0 or lz4x:1 for runtime tasks, defl:3 elsewhere.

# Some metrics (YMMV)

```
C:\prj\stdarc>test.exe --benchmark --0 \prj\enwik8 | sort
Y 100000000 ->    63186893 defl.0 63.19% c:1.670s d:0.418s \prj\enwik8
Y 100000000 ->    58555380  ppp.0 58.56% c:0.684s d:0.393s \prj\enwik8
Y 100000000 ->    56013654 lzp1.0 56.01% c:0.778s d:0.508s \prj\enwik8
Y 100000000 ->    52650827 lz4x.0 52.65% c:0.688s d:0.097s \prj\enwik8
Y 100000000 ->    51923215  ulz.0 51.92% c:0.602s d:0.089s \prj\enwik8
Y 100000000 ->    49061221 lzss.0 49.06% c:13.96s d:0.344s \prj\enwik8
Y 100000000 ->    47874204 lzw3.0 47.87% c:1.069s d:0.286s \prj\enwik8
Y 100000000 ->    47620393 crsh.0 47.62% c:1.797s d:0.449s \prj\enwik8
Y 100000000 ->    39963844 lzma.0 39.96% c:13.53s d:2.680s \prj\enwik8
Y 100000000 ->    30056093 balz.0 30.06% c:10.19s d:5.109s \prj\enwik8
Y 100000000 ->    27700930  bcm.0 27.70% c:13.59s d:11.71s \prj\enwik8
80.64s
```

```
C:\prj\stdarc>test.exe --benchmark --9 \prj\enwik8 | sort
Y 100000000 ->    58555380  ppp.9 58.56% c:0.688s d:0.394s \prj\enwik8
Y 100000000 ->    56013654 lzp1.9 56.01% c:0.769s d:0.507s \prj\enwik8
Y 100000000 ->    49061221 lzss.9 49.06% c:13.58s d:0.334s \prj\enwik8
Y 100000000 ->    47874204 lzw3.9 47.87% c:1.064s d:0.286s \prj\enwik8
Y 100000000 ->    46021464 lz4x.9 46.02% c:1.293s d:0.095s \prj\enwik8
Y 100000000 ->    40253295  ulz.9 40.25% c:11.78s d:0.097s \prj\enwik8
Y 100000000 ->    36460096 defl.9 36.46% c:6.372s d:0.387s \prj\enwik8
Y 100000000 ->    32038271 crsh.9 32.04% c:107.9s d:0.343s \prj\enwik8
Y 100000000 ->    28232820 balz.9 28.23% c:32.63s d:4.855s \prj\enwik8
Y 100000000 ->    25801721 lzma.9 25.80% c:93.58s d:1.413s \prj\enwik8
Y 100000000 ->    20789660  bcm.9 20.79% c:19.72s d:17.84s \prj\enwik8
315.97s
```

# Project goals
- C.
- Drop & use.
- Small codebase.
- Good enough performance.
- De/compression for memory blocks, files and archives.
- Public domain or unlicensed source code. No attribution required.

# Todo
- Test other environments. Currently VS2019+GCC+Clang only.
- Add small compression filters (like E8E9).
- ~~Add vfs and archive mounting support.~~
- ~~Add small file archivers (pak, zip, tar).~~
- Optimize LZSS. Could compressor be nearly as performant as LZJB with some effort?
- Optimize BALZ. Could those in-mem bufferings be skipped?
- Consider ZPAQ (large C++), ~~BCM~~, ~~yalz77~~. ~~LZJB maybe?~~
- Bugfixes and cleaning up.

# Links
- Ilya Muravyov: bcm, balz, crush, ulz, lz4x - https://github.com/encode84
- Rich Geldreich: miniz - https://github.com/richgel999
- Igor Pavlov: lzma - https://www.7-zip.org/
- Charles Bloom: lzp1 - https://www.cbloom.com/src/index_lz.html
- Ross Williams: lzrw3a - http://ross.net/compression/lzrw3a.html
- Haruhiko Okumura: lzss - https://oku.edu.mie-u.ac.jp/~okumura/compression/
- Dave Rand: ppp - https://tools.ietf.org/html/rfc1978
- Joonas Pihlajamaa: junzip - https://github.com/jokkebk/JUnzip/
