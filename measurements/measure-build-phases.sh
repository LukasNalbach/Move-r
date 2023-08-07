while getopts p:t: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
        t) t=${OPTARG};;
    esac
done

for ((i_p=1; i_p<=$p; i_p*=2)); do
    ../build/cli/move-r-build -o "indexes/$t-8" -p $i_p -a 8 -m_idx "results/results-move-r-build-phases.txt" -m_mds "results/results-mds-build-phases.txt" "texts/$t"
done