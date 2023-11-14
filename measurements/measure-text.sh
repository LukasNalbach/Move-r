while getopts p:t: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
        t) t=${OPTARG};;
    esac
done

cd ..

./build/cli/move-r-bench -m measurements/results/results.txt measurements/texts/$t measurements/patterns/$t-patterns-bal measurements/patterns/$t-patterns-phi $p
./build/cli/move-r-bench -a -m measurements/results/results-a.txt measurements/texts/$t measurements/patterns/$t-patterns-bal measurements/patterns/$t-patterns-phi $p

cd measurements

for ((i_p=1; i_p<=$p; i_p*=2)); do
    ../build/cli/move-r-build -o "indexes/$t-8" -p $i_p -a 8 -m_idx "results/results-move-r-build-phases.txt" -m_mds "results/results-mds-build-phases.txt" "texts/$t"
done