while getopts p:t: flag
do
    case "${flag}" in
        p) p=${OPTARG};;
        t) t=${OPTARG};;
    esac
done

./measure-build-phases.sh -p $p -t $t
./measure-a.sh -p $p -t $t