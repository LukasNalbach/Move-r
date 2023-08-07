while getopts p:t: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
        t) t=${OPTARG};;
    esac
done

cd ..

./build/cli/move-r-bench -a -m measurements/results/results-a.txt measurements/texts/$t measurements/patterns/$t-patterns-bal measurements/patterns/$t-patterns-phi $p

cd measurements