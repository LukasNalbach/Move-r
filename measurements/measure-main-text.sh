while getopts p:t: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
        t) t=${OPTARG};;
    esac
done

cd ..

./build/cli/move-r-bench -m measurements/results/results-main.txt measurements/texts/$t measurements/patterns/$t-patterns-bal measurements/patterns/$t-patterns-phi $p
./build/cli/move-r-bench -sa -m measurements/results/results-main.txt measurements/texts/$t $p

cd measurements