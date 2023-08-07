while getopts p: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
    esac
done

./measure-main.sh -p $p
./measure-appendix.sh -p $p