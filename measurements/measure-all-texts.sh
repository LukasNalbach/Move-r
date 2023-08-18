while getopts p: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
    esac
done

./measure-text.sh -p $p -t einstein.de.txt
./measure-text.sh -p $p -t Horspool
./measure-text.sh -p $p -t dewiki.8GiB
./measure-text.sh -p $p -t dna.001.1
./measure-text.sh -p $p -t english.200MB
./measure-text.sh -p $p -t dna.200MB