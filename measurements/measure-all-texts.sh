while getopts p: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
    esac
done

./measure-text.sh -p $p -t einstein.en.txt
./measure-text.sh -p $p -t english
./measure-text.sh -p $p -t dewiki.64GiB
./measure-text.sh -p $p -t chr.19.1000
./measure-text.sh -p $p -t sars.cov.2