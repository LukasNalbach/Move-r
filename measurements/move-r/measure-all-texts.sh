while getopts p: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
    esac
done

./measure-text.sh -p $p -t einstein.en.txt
./measure-text.sh -p $p -t english
./measure-text.sh -p $p -t dewiki
./measure-text.sh -p $p -t sars2
./measure-text.sh -p $p -t chr19