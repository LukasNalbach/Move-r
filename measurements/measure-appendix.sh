while getopts p: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
    esac
done

./measure-appendix-text.sh -p $p -t world_leaders
./measure-appendix-text.sh -p $p -t EDGAR
./measure-appendix-text.sh -p $p -t einstein.de.txt
./measure-appendix-text.sh -p $p -t Horspool
./measure-appendix-text.sh -p $p -t dewiki.8GiB
./measure-appendix-text.sh -p $p -t dna.200MB
./measure-appendix-text.sh -p $p -t Quicksort
./measure-appendix-text.sh -p $p -t english