# Approximate pattern matching: move_rb against sassy, on chr19 and sars2

Measured 2026-09-10/11 on the 9800X3D (8 cores / 16 threads, 96 MiB L3, 62 GiB RAM, WSL2).
All queries single-threaded.

`move_rb` is an index. [sassy](https://github.com/RagnarGrootKoerkamp/sassy) 0.2.6 (Beeloo & Groot
Koerkamp) has no index: it scans the whole text for each query with a SIMD-tiled bitparallel Myers
algorithm. Both search the same text with the same patterns.

## The two texts

Both texts are FASTA files and were indexed in FASTA mode (`move-rb build -f`): the headers and line
breaks are removed and each run of non-`ACGT` bases is replaced by one separator symbol.
`fasta-to-text` wrote the resulting text to a file. sassy searches this file and the patterns were
sampled from it. Its size is equal to `n` of the index.

| | sars2 | chr19 |
|---|---|---|
| FASTA | 88.20 GB, 2 904 807 sequences | 60.11 GB, 1000 sequences |
| text after FASTA mode (`n`) | 84 593 139 404 | 55 805 122 010 |
| alphabet | 6 (`ACGT`, separator, terminator) | 6 |
| runs `r` (forward / backward) | 84 141 607 / 83 924 917 | 45 904 750 / 45 918 212 |
| `n/r` | 1005 / 1007 | 1215 |
| index (`locate_move`) | **6113 MB** (7.3 % of the text) | **3093 MB** (5.5 % of the text) |
| construction | 24:30 min, peak 7304 MB | 15:58 min, peak 3649 MB |

Both were built with Big-BWT (`-c bigbwt`), 16 threads, `a = 8`, `-s locate_move`.

## Differences between move_rb and sassy

1. **sassy supports only locate queries under the edit distance.** The tables for count and for the
   Hamming distance below contain only move_rb.
2. **Its `dna` profile is not exact on these texts.** It encodes a byte `c` as `(c >> 1) & 3`, i.e.,
   every byte is mapped to a base; the FASTA separator `0x7F` is mapped to `G`. The measurements
   below use the `ascii` profile, which is exact. `dna` is 7-20 % faster and reports about 0.1 %
   additional matches that span a separator (measured on sars2 with one pattern per cell):

   | k, m | `dna` | `ascii` | occurrences `dna` / `ascii` |
   |---|---|---|---|
   | 1, 40 | 17 s | 21 s | 2 757 914 / 2 755 209 |
   | 1, 160 | 20 s | 25 s | 2 805 470 / 2 805 074 |
   | 1, 640 | 51 s | 56 s | 2 193 706 / 2 191 183 |
   | 7, 40 | 28 s | 33 s | 2 885 300 / 2 885 001 |
   | 7, 160 | 30 s | 35 s | 2 833 077 / 2 832 587 |
   | 7, 640 | 1:07 min | 1:12 min | 2 762 543 / 2 762 196 |

3. **The two report different sets of occurrences.** sassy's `search` reports the rightmost position
   of each local minimum of the end position; move_rb reports the occurrences by their starting
   positions. If a pattern has few occurrences, both report the same number (see the chr19 table:
   1000 and 1000, 795 and 795, 736 and 736); if it has many, move_rb reports 2-20x more. The tables
   below therefore contain the number of occurrences per pattern and the time per reported
   occurrence in addition to the time per query.

## locate, edit distance: move_rb and sassy

Both ran the same patterns once: the first two patterns of each file
(`--max-patterns 2 --time 0` for move_rb, `--sassy-patterns 2`). The patterns of a file differ
widely in their numbers of occurrences, so on chr19 the average over the first two patterns differs
from the average over all patterns.

`/occ` is the time per reported occurrence. `sassy/move_rb` is the time of sassy divided by the time
of move_rb.

### sars2 (n = 84.6e9)

```
--- locate / edit ------------------------------------------------------------
  k     m | move_rb_move      occ/pat     /occ | sassy_ascii      occ/pat     /occ | sassy/move_rb
  1    10 |      1.43 s   11,368,444 126.0 ns |     21.48 s    8,500,092   2.5 us |          15.0
  1    20 |    695.0 ms    5,644,052 123.1 ns |     21.55 s    2,826,920   7.6 us |          31.0
  1    40 |    692.6 ms    5,600,951 123.7 ns |     21.94 s    2,805,179   7.8 us |          31.7
  1    80 |    686.8 ms    5,637,681 121.8 ns |     23.02 s    2,848,961   8.1 us |          33.5
  1   160 |    663.2 ms    5,434,802 122.0 ns |     25.48 s    2,831,596   9.0 us |          38.4
  1   320 |    400.4 ms    3,149,798 127.1 ns |     30.19 s    1,696,434  17.8 us |          75.4
  1   640 |    437.3 ms    3,382,315 129.3 ns |     57.64 s    2,312,128  24.9 us |         131.8
  1  1280 |    286.5 ms    1,871,919 153.1 ns |     94.37 s    1,194,742  79.0 us |         329.4
  2    20 |      1.07 s    8,607,840 124.3 ns |     27.12 s    2,877,186   9.4 us |          25.4
  2    40 |    868.6 ms    7,144,619 121.6 ns |     26.02 s    2,861,034   9.1 us |          30.0
  2    80 |    898.8 ms    7,426,836 121.0 ns |     28.42 s    2,825,644  10.1 us |          31.6
  2   160 |    655.9 ms    5,106,715 128.4 ns |     30.25 s    2,732,933  11.1 us |          46.1
  2   320 |    771.7 ms    6,090,062 126.7 ns |     37.72 s    2,316,115  16.3 us |          48.9
  2   640 |    509.4 ms    3,594,128 141.7 ns |     52.82 s    1,559,020  33.9 us |         103.7
  2  1280 |    650.3 ms    3,188,089 204.0 ns |    129.56 s    1,690,981  76.6 us |         199.2
  4    20 |      1.72 s   13,955,176 123.2 ns |     27.81 s    2,831,610   9.8 us |          16.2
  4    40 |      1.16 s    9,339,024 124.5 ns |     29.00 s    2,696,354  10.8 us |          24.9
  4    80 |      1.56 s   12,663,155 123.3 ns |     30.17 s    2,838,743  10.6 us |          19.3
  4   160 |      1.39 s   11,159,852 125.0 ns |     32.64 s    2,763,082  11.8 us |          23.4
  4   320 |      1.39 s   10,594,232 131.2 ns |     41.22 s    2,806,490  14.7 us |          29.7
  4   640 |      1.06 s    6,505,010 163.0 ns |     66.22 s    2,389,984  27.7 us |          62.5
  4  1280 |      1.39 s    4,723,757 293.2 ns |    153.15 s    2,159,605  70.9 us |         110.6
  7    40 |      2.50 s   19,663,162 126.9 ns |     33.15 s    2,850,779  11.6 us |          13.3
  7    80 |      2.75 s   21,410,732 128.4 ns |     33.46 s    2,809,621  11.9 us |          12.2
  7   160 |      2.15 s   16,976,025 126.7 ns |     36.12 s    2,835,792  12.7 us |          16.8
  7   320 |      1.67 s   11,657,288 143.6 ns |     45.35 s    2,760,509  16.4 us |          27.1
  7   640 |      2.52 s   12,897,534 195.8 ns |     72.13 s    2,739,010  26.3 us |          28.6
  7  1280 |      4.96 s    8,929,710 555.1 ns |    163.40 s    2,619,410  62.4 us |          33.0
 10    40 |      2.85 s   22,211,247 128.2 ns |     37.37 s    2,795,988  13.4 us |          13.1
 10    80 |      2.50 s   18,571,737 134.8 ns |     37.87 s    2,775,457  13.6 us |          15.1
 10   160 |      2.17 s   16,011,878 135.4 ns |     40.23 s    2,636,038  15.3 us |          18.6
 10   320 |      2.34 s   17,240,678 135.9 ns |     48.92 s    2,801,970  17.5 us |          20.9
 10   640 |      3.52 s   19,364,125 181.7 ns |     75.25 s    2,744,610  27.4 us |          21.4
 10  1280 |      6.39 s   10,152,579 629.6 ns |    143.68 s    2,037,962  70.5 us |          22.5
 13    40 |      5.75 s   39,254,407 146.6 ns |     40.72 s    8,465,492   4.8 us |           7.1
 13    80 |      1.97 s   15,412,432 127.8 ns |     41.64 s    2,879,142  14.5 us |          21.1
 13   160 |      2.35 s   17,682,202 132.9 ns |     44.80 s    2,831,611  15.8 us |          19.1
 13   320 |      2.75 s   18,173,572 151.1 ns |     51.89 s    2,757,546  18.8 us |          18.9
 13   640 |      3.01 s   10,419,122 289.0 ns |     73.65 s    2,150,806  34.2 us |          24.5
 13  1280 |      6.87 s    9,634,984 712.9 ns |    130.59 s    1,732,014  75.4 us |          19.0
```

### chr19 (n = 55.8e9)

```
--- locate / edit ------------------------------------------------------------
  k     m | move_rb_move      occ/pat     /occ | sassy_ascii      occ/pat     /occ | sassy/move_rb
  1    10 |      1.79 s   19,668,434  91.1 ns |     16.25 s   16,050,158   1.0 us |           9.1
  1    20 |    405.0 us        1,000 405.0 ns |     13.62 s        1,000  13.6 ms |         33.6k
  1    40 |    309.1 us        1,500 206.0 ns |     13.56 s        1,000  13.6 ms |         43.9k
  1    80 |    314.5 us        1,000 314.5 ns |     13.71 s        1,000  13.7 ms |         43.6k
  1   160 |    321.7 us        1,000 321.7 ns |     13.51 s        1,000  13.5 ms |         42.0k
  1   320 |    441.8 us          795 555.7 ns |     14.59 s          795  18.4 ms |         33.0k
  1   640 |    658.6 us          808 815.2 ns |     13.57 s          808  16.8 ms |         20.6k
  1  1280 |      1.2 ms          736   1.6 us |     13.74 s          736  18.7 ms |         11.6k
  2    10 |      3.82 s   38,769,983  98.7 ns |     18.90 s   29,310,009 644.7 ns |           4.9
  2    20 |      1.20 s   12,051,788  99.3 ns |     18.64 s    7,182,612   2.6 us |          15.6
  2    40 |    396.2 us        1,000 396.2 ns |     16.37 s        1,000  16.4 ms |         41.3k
  2    80 |    336.3 us        1,000 336.3 ns |     17.72 s        1,000  17.7 ms |         52.7k
  2   160 |    394.4 us        1,000 394.4 ns |     17.86 s        1,000  17.9 ms |         45.3k
  2   320 |    475.6 us          802 593.1 ns |     16.94 s          802  21.1 ms |         35.6k
  2   640 |    670.9 us          668   1.0 us |     16.77 s          668  25.1 ms |         25.0k
  2  1280 |      1.4 ms          637   2.2 us |     16.68 s          637  26.2 ms |         11.7k
  4    20 |    135.1 ms      722,880 187.0 ns |     18.49 s      597,428  30.9 us |         136.8
  4    40 |    789.4 us        1,001 788.6 ns |     18.78 s        1,000  18.8 ms |         23.8k
  4    80 |      1.6 ms        4,954 321.5 ns |     19.68 s        4,454   4.4 ms |         12.4k
  4   160 |    389.1 us        1,000 389.1 ns |     18.95 s        1,000  18.9 ms |         48.7k
  4   320 |    487.6 us        1,000 487.6 ns |     18.75 s        1,000  18.8 ms |         38.5k
  4   640 |      1.4 ms        3,492 401.9 ns |     18.94 s        3,492   5.4 ms |         13.5k
  4  1280 |      2.8 ms        1,000   2.8 us |     18.83 s        1,000  18.8 ms |          6.7k
  7    40 |     10.7 ms        2,000   5.3 us |     21.37 s        1,000  21.4 ms |          2.0k
  7    80 |      8.3 ms       75,038 110.8 ns |     21.50 s       36,920 582.4 us |          2.6k
  7   160 |    497.0 us        2,000 248.5 ns |     21.28 s        2,000  10.6 ms |         42.8k
  7   320 |    610.1 us        1,000 610.1 ns |     22.78 s        1,000  22.8 ms |         37.3k
  7   640 |    892.0 us        1,000 892.0 ns |     23.30 s        1,000  23.3 ms |         26.1k
  7  1280 |      3.6 ms        1,000   3.6 us |     21.23 s        1,000  21.2 ms |          5.9k
 10    40 |      1.06 s        1,535 690.5 us |     23.95 s        1,501  16.0 ms |          22.6
 10    80 |    411.5 ms    1,768,474 232.7 ns |     28.40 s    1,041,950  27.3 us |          69.0
 10   160 |    836.9 us        1,000 836.9 ns |     23.99 s        1,228  19.5 ms |         28.7k
 10   320 |    756.4 us        1,000 756.4 ns |     24.04 s        1,000  24.0 ms |         31.8k
 10   640 |      1.8 ms        1,000   1.8 us |     26.10 s        1,000  26.1 ms |         14.5k
 10  1280 |      4.8 ms          928   5.2 us |     24.19 s          928  26.0 ms |          5.0k
 13    40 |      7.94 s      676,586  11.7 us |     26.29 s      567,065  46.4 us |           3.3
 13    80 |    987.6 ms       12,206  80.9 us |     26.67 s        6,320   4.2 ms |          27.0
 13   160 |    748.4 us        1,000 748.4 ns |     26.88 s        1,000  26.9 ms |         35.9k
 13   320 |    916.8 us        1,003 914.0 ns |     26.61 s        1,000  26.6 ms |         29.0k
 13   640 |      2.5 ms        1,000   2.5 us |     28.83 s        1,000  28.8 ms |         11.4k
 13  1280 |      7.5 ms        1,000   7.5 us |     32.10 s        1,000  32.1 ms |          4.3k
```

**Summary.** move_rb takes **~120-150 ns per reported occurrence on sars2** and 0.1-1 us on chr19,
almost independent of k and m; its time depends on the number of occurrences. The time of sassy
depends on the length of the text: 13.5 to 163 s per query, independent of the number of
occurrences, i.e., 0.6 us to 32 ms per occurrence. The ratio is smallest where the number of
occurrences is largest: 7.1x on sars2 at k=13, m=40 (39 million occurrences per pattern), and
43 900x on chr19 at k=1, m=40 (1500 occurrences per pattern).

## Scan rate of sassy

`time_search` is the time spent in sassy, without reading the text. The next 256 MB chunk is read
by a second thread while the current chunk is searched; the wall time and the search time differ by
less than 0.5 % in all cells. (Read times measured on this machine are not meaningful: WSL2 serves
O_DIRECT reads from the cache of the Windows host.)

```
sars2:
    m k=1       k=2       k=4       k=7       k=10      k=13     
   10 3.95      -         -         -         -         -           [GB/s]
   20 3.94      3.13      3.05      -         -         -           [GB/s]
   40 3.87      3.26      2.92      2.56      2.27      2.08        [GB/s]
   80 3.69      2.98      2.81      2.53      2.24      2.03        [GB/s]
  160 3.33      2.80      2.60      2.35      2.11      1.89        [GB/s]
  320 2.81      2.25      2.06      1.87      1.73      1.63        [GB/s]
  640 1.47      1.60      1.28      1.17      1.13      1.15        [GB/s]
 1280 0.90      0.65      0.55      0.52      0.59      0.65        [GB/s]
chr19:
    m k=1       k=2       k=4       k=7       k=10      k=13     
   10 3.45      2.96      -         -         -         -           [GB/s]
   20 4.12      3.00      3.03      -         -         -           [GB/s]
   40 4.13      3.42      2.98      2.62      2.34      2.13        [GB/s]
   80 4.09      3.16      2.85      2.60      1.97      2.10        [GB/s]
  160 4.15      3.14      2.96      2.63      2.33      2.08        [GB/s]
  320 3.84      3.31      2.99      2.46      2.33      2.10        [GB/s]
  640 4.13      3.34      2.96      2.40      2.14      1.94        [GB/s]
 1280 4.08      3.36      2.97      2.64      2.31      1.74        [GB/s]
```

**On chr19 the rate does not depend on m**: 4.1 GB/s at m = 10 and 4.1 GB/s at m = 1280. The
pattern matches at almost no position, so the computation of a bitparallel column stops after its
first words. **On sars2 the rate decreases with m**, from 3.9 GB/s to 0.90 GB/s. sars2 consists of
2.9 million nearly identical genomes, so the pattern nearly matches at every position, the whole
column is computed and the time is O(n * m / 64).

## count and locate under Hamming distance (move_rb only)

sassy supports neither operation. These were measured with the whole pattern sets; each file was
replayed for at least 10 s.

### sars2

```
--- count / hamming ------------------------------------------------------------
  k     m | move_rb_move      occ/pat     /occ
  1    10 |      7.5 us    8,332,817   0.0 ns
  1    20 |     33.7 us    2,816,921   0.0 ns
  1    40 |     97.9 us    2,785,145   0.0 ns
  1    80 |    284.6 us    2,722,416   0.1 ns
  1   160 |    869.3 us    2,544,931   0.3 ns
  1   320 |      2.7 ms    2,341,052   1.2 ns
  1   640 |      8.5 ms    1,837,101   4.6 ns
  1  1280 |     26.0 ms    1,322,988  19.7 ns
  2    10 |     84.2 us   67,676,669   0.0 ns
  2    20 |    112.9 us    2,825,846   0.0 ns
  2    40 |    243.4 us    2,795,998   0.1 ns
  2    80 |    838.2 us    2,751,511   0.3 ns
  2   160 |      3.2 ms    2,671,390   1.2 ns
  2   320 |     12.5 ms    2,422,987   5.2 ns
  2   640 |     40.7 ms    1,784,364  22.8 ns
  2  1280 |    149.4 ms    1,828,843  81.7 ns
  4    10 |      2.6 ms 2,343,797,095   0.0 ns
  4    20 |    604.0 us    2,949,711   0.2 ns
  4    40 |    548.7 us    2,810,475   0.2 ns
  4    80 |      2.0 ms    2,782,107   0.7 ns
  4   160 |      8.7 ms    2,772,327   3.1 ns
  4   320 |     35.8 ms    2,781,387  12.9 ns
  4   640 |    194.6 ms    2,616,610  74.4 ns
  4  1280 |    934.0 ms    2,497,678 373.9 ns
  7    10 |     74.3 ms 42,313,597,567   0.0 ns
  7    20 |     12.7 ms   36,585,819   0.3 ns
  7    40 |    932.9 us    2,809,415   0.3 ns
  7    80 |      3.1 ms    2,760,497   1.1 ns
  7   160 |     16.4 ms    2,763,414   5.9 ns
  7   320 |     68.2 ms    2,447,449  27.8 ns
  7   640 |    406.1 ms    2,600,806 156.1 ns
  7  1280 |      2.34 s    2,587,555 904.0 ns
 10    20 |    456.9 ms 2,212,668,184   0.2 ns
 10    40 |      4.7 ms    2,823,358   1.7 ns
 10    80 |      4.4 ms    2,798,772   1.6 ns
 10   160 |     20.1 ms    2,683,731   7.5 ns
 10   320 |     90.2 ms    2,791,520  32.3 ns
 10   640 |    433.9 ms    2,715,392 159.8 ns
 10  1280 |      2.38 s    2,588,917 919.5 ns
 13    20 |      1.63 s 21,907,430,917   0.1 ns
 13    40 |     23.4 ms    2,789,612   8.4 ns
 13    80 |      5.5 ms    2,785,916   2.0 ns
 13   160 |     24.3 ms    2,767,022   8.8 ns
 13   320 |    110.7 ms    2,677,621  41.4 ns
 13   640 |    544.2 ms    2,734,854 199.0 ns
 13  1280 |      3.86 s    1,689,424   2.3 us
```

```
--- locate / hamming ------------------------------------------------------------
  k     m | move_rb_move      occ/pat     /occ
  1    10 |    461.1 ms    4,253,992 108.4 ns
  1    20 |    302.7 ms    2,831,572 106.9 ns
  1    40 |    312.1 ms    2,868,580 108.8 ns
  1    80 |    300.2 ms    2,756,311 108.9 ns
  1   160 |    304.5 ms    2,792,335 109.1 ns
  1   320 |    168.9 ms    1,531,738 110.2 ns
  1   640 |    225.1 ms    2,018,651 111.5 ns
  1  1280 |    217.4 ms    1,709,902 127.1 ns
  2    10 |      5.14 s   46,205,452 111.2 ns
  2    20 |    302.1 ms    2,853,930 105.9 ns
  2    40 |    298.3 ms    2,816,957 105.9 ns
  2    80 |    305.3 ms    2,828,776 107.9 ns
  2   160 |    309.9 ms    2,846,877 108.9 ns
  2   320 |    310.1 ms    2,697,134 115.0 ns
  2   640 |    309.7 ms    2,415,029 128.3 ns
  2  1280 |    447.1 ms    2,263,436 197.5 ns
  4    20 |    442.2 ms    4,148,394 106.6 ns
  4    40 |    304.5 ms    2,875,961 105.9 ns
  4    80 |    299.9 ms    2,795,005 107.3 ns
  4   160 |    305.4 ms    2,767,553 110.4 ns
  4   320 |    292.3 ms    2,413,640 121.1 ns
  4   640 |    392.8 ms    2,004,945 195.9 ns
  4  1280 |      1.09 s    2,001,292 544.8 ns
  7    20 |      3.54 s   32,217,568 109.8 ns
  7    40 |    297.2 ms    2,778,799 106.9 ns
  7    80 |    297.8 ms    2,788,155 106.8 ns
  7   160 |    325.1 ms    2,775,229 117.1 ns
  7   320 |    453.0 ms    2,747,147 164.9 ns
  7   640 |    556.3 ms    2,735,250 203.4 ns
  7  1280 |      1.80 s    1,550,368   1.2 us
 10    40 |    298.2 ms    2,716,933 109.7 ns
 10    80 |    309.2 ms    2,847,595 108.6 ns
 10   160 |    330.5 ms    2,738,436 120.7 ns
 10   320 |    389.0 ms    2,665,526 146.0 ns
 10   640 |    736.5 ms    2,226,478 330.8 ns
 10  1280 |      3.00 s    2,641,418   1.1 us
 13    40 |    335.2 ms    2,866,130 116.9 ns
 13    80 |    311.5 ms    2,849,425 109.3 ns
 13   160 |    315.4 ms    2,659,721 118.6 ns
 13   320 |    488.5 ms    2,726,594 179.2 ns
 13   640 |    849.0 ms    2,709,160 313.4 ns
 13  1280 |      4.28 s    2,664,534   1.6 us
```

### chr19

```
--- count / hamming ------------------------------------------------------------
  k     m | move_rb_move      occ/pat     /occ
  1    10 |      7.4 us    7,381,107   0.0 ns
  1    20 |     25.6 us      506,943   0.1 ns
  1    40 |     22.7 us       15,877   1.4 ns
  1    80 |     30.4 us        1,131  26.9 ns
  1   160 |     57.8 us          994  58.1 ns
  1   320 |    115.1 us          918 125.5 ns
  1   640 |    240.4 us          794 302.7 ns
  1  1280 |    501.3 us          601 833.8 ns
  2    10 |     94.3 us   45,794,566   0.0 ns
  2    20 |    281.3 us    1,266,429   0.2 ns
  2    40 |     90.0 us       60,452   1.5 ns
  2    80 |     53.4 us        1,371  38.9 ns
  2   160 |     80.2 us        1,026  78.2 ns
  2   320 |    172.1 us          962 178.8 ns
  2   640 |    396.5 us          872 454.6 ns
  2  1280 |    915.7 us          717   1.3 us
  4    10 |      3.7 ms 1,327,235,103   0.0 ns
  4    20 |      5.7 ms    3,979,672   1.4 ns
  4    40 |      1.0 ms      369,184   2.7 ns
  4    80 |    200.0 us        3,871  51.7 ns
  4   160 |    132.8 us        1,075 123.6 ns
  4   320 |    283.7 us          983 288.7 ns
  4   640 |    725.5 us          920 789.0 ns
  4  1280 |      1.9 ms          782   2.5 us
  7    10 |    141.8 ms 26,239,916,923   0.0 ns
  7    20 |    219.7 ms   26,355,517   8.3 ns
  7    40 |      8.0 ms    1,699,682   4.7 ns
  7    80 |      1.5 ms       24,168  60.5 ns
  7   160 |    274.7 us        1,279 214.8 ns
  7   320 |    395.0 us          985 401.1 ns
  7   640 |      1.1 ms          929   1.2 us
  7  1280 |      3.2 ms          828   3.9 us
 10    20 |      8.98 s  840,080,267  10.7 ns
 10    40 |     86.9 ms    1,008,686  86.2 ns
 10    80 |      8.1 ms      152,742  52.7 ns
 10   160 |    972.2 us        1,795 541.6 ns
 10   320 |    539.0 us        1,045 515.6 ns
 10   640 |      1.4 ms          953   1.4 us
 10  1280 |      4.0 ms          834   4.8 us
 13    40 |    493.3 ms    1,042,190 473.3 ns
 13    80 |     12.5 ms       96,750 129.7 ns
 13   160 |      3.3 ms        3,367 992.3 ns
 13   320 |    779.5 us        1,059 736.3 ns
 13   640 |      1.7 ms          943   1.8 us
 13  1280 |      5.3 ms          834   6.3 us
```

```
--- locate / hamming ------------------------------------------------------------
  k     m | move_rb_move      occ/pat     /occ
  1    10 |    179.6 ms    2,504,346  71.7 ns
  1    20 |     57.9 ms      853,074  67.8 ns
  1    40 |      1.1 ms       13,403  80.0 ns
  1    80 |    137.5 us        1,139 120.7 ns
  1   160 |    157.3 us        1,009 155.9 ns
  1   320 |    218.7 us          916 238.8 ns
  1   640 |    345.4 us          793 435.4 ns
  1  1280 |    604.8 us          627 964.6 ns
  2    10 |      4.73 s   68,036,964  69.4 ns
  2    20 |    174.1 ms    2,234,420  77.9 ns
  2    40 |      4.0 ms       48,438  82.6 ns
  2    80 |    191.4 us        1,426 134.2 ns
  2   160 |    190.0 us        1,028 184.7 ns
  2   320 |    286.4 us          967 296.1 ns
  2   640 |    520.4 us          860 605.1 ns
  2  1280 |      1.1 ms          726   1.5 us
  4    20 |    682.1 ms    8,584,249  79.5 ns
  4    40 |    115.2 ms    1,381,136  83.4 ns
  4    80 |    890.2 us        6,735 132.2 ns
  4   160 |    256.6 us        1,105 232.4 ns
  4   320 |    418.1 us          969 431.7 ns
  4   640 |    920.4 us          916   1.0 us
  4  1280 |      2.2 ms          799   2.8 us
  7    20 |      1.45 s   13,782,341 104.9 ns
  7    40 |     91.1 ms    1,091,813  83.5 ns
  7    80 |      7.4 ms       52,577 141.7 ns
  7   160 |    432.4 us        1,277 338.6 ns
  7   320 |    540.7 us          993 544.5 ns
  7   640 |      1.3 ms          925   1.4 us
  7  1280 |      3.7 ms          872   4.2 us
 10    40 |    341.6 ms    2,998,436 113.9 ns
 10    80 |     28.9 ms      200,898 143.8 ns
 10   160 |      1.4 ms        2,392 602.7 ns
 10   320 |    690.1 us        1,010 683.0 ns
 10   640 |      1.6 ms          935   1.7 us
 10  1280 |      4.5 ms          851   5.3 us
 13    40 |    446.6 ms      360,936   1.2 us
 13    80 |     61.6 ms      389,192 158.3 ns
 13   160 |      4.2 ms        4,107   1.0 us
 13   320 |    910.9 us        1,022 891.3 ns
 13   640 |      1.9 ms          951   2.0 us
 13  1280 |      5.7 ms          846   6.7 us
```

## Reproducing

```
measurements/move-rb/measure-apm-sassy.sh -f ~/texts/sars2 -S 2 -P ascii -R 10
measurements/move-rb/measure-apm-sassy.sh -f ~/texts/chr19 -S 2 -P ascii -R 10
```

The script builds the index in FASTA mode, writes the plain text with `fasta-to-text`, generates the
pattern sets with `move-rb-gen-apm-queries` (k = 1,2,4,7,10,13, m = 10,20,40,80,160,320,640,1280,
seed 42, patterns without the separator) and runs `move-rb-bench-apm`. Of the 48 (k, m) cells, the
generator omits those in which no query finishes within its timeout of 10 s: 2-3 for count/Hamming,
6 for locate/Hamming and 8 (sars2) resp. 7 (chr19) for locate/edit. All of them have a high error
rate (k >= 2 at m = 10, k >= 7 at m = 20). 40 resp. 41 locate/edit cells remain for the comparison.

The tables comparing move_rb and sassy were measured in a second run with the same patterns:

```
build/bench/move-rb-bench-apm --metric edit --only locate --time 0 \
    --max-patterns 2 \
    --move-rb ~/texts/<text>.move-rb-locate \
    -m results-apm-head2head.txt \
    <text> ~/texts/<text>.patterns-apm/
```

`RESULT` lines: `results-apm-sassy.txt` (whole pattern sets, both texts, with sassy),
`results-apm-head2head.txt` (move_rb on the first two patterns) and
`results-apm-sassy-profiles.txt` (comparison of the `dna` and `ascii` profiles).
