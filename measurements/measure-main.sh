while getopts p: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
    esac
done

./measure-main-text.sh -p $p -t einstein.de.txt
./measure-main-text.sh -p $p -t Horspool
./measure-main-text.sh -p $p -t dewiki.8GiB
./measure-main-text.sh -p $p -t dna.200MB