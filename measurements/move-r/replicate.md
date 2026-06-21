#### How to replicate the measurements
1. Build the project with `MOVE_R_BUILD_BENCH_CLI` set to `ON`.
2. Download and decompress the texts:
- [einstein.en.txt](https://pizzachili.dcc.uchile.cl/repcorpus/real/einstein.en.txt.7z)
- [english](https://pizzachili.dcc.uchile.cl/texts/nlang/english.gz)
- [chr19](https://drive.google.com/file/d/1GrCHHcc3zH56Q-c6WbI1N6qOh0sBD5DO/view?usp=sharing)
- [dewiki](https://drive.google.com/file/d/1GqvkN0FH6dkSxHZCXFPOr7I1iOBenUIZ/view?usp=sharing)
- [sars2](https://drive.google.com/file/d/134fLOpY1_3dFTdSSc_vW2qai4yKuyl2W/view?usp=sharing)
3. Place the texts into the folder `measurements/texts/`.
4. Navigate to the folder `measurements/`.
5. Run the measurements:
To measure all texts, run `./measure-all-texts.sh -p "num_threads"`.
To measure a single text, run `./measure-text.sh -t "text_name" -p "num_threads"`.
`"num_threads"` is the maximum number of threads to use.
6. To measure other texts, generate two sets of patterns with `move-r-patterns` using
the same naming scheme and place them into the folder `measurements/patterns/`.
7. The results are written to files in the folder `measurements/results/`. To import
them into LaTeX, use [sqlplot-tools](https://github.com/bingmann/sqlplot-tools).