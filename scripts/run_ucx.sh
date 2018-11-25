#!/bin/bash

OMPI_HOME=/labhome/xinz/workplace/ompi/install
UCX_HOME=/labhome/xinz/workplace/ucx/install

export PATH=$OMPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$UCX_HOME/lib:$LD_LIBRARY_PATH

MSG_RANGE=1:8192:2
WIN_SIZE=512
PE=2
BLOCKING=
BINDING=
LOOP=100
WLOOP=10
TYPE=p
UCX_DEV=mlx5_0:1
UCX_TRANS=rc_x
VERBOSE=0

function usage()
{
    echo -e "\t-h\t--help"
    echo -e "\t-m\tmessage range {lb:hb:stride}"
    echo -e "\t-w\twindow size (default:512)"
    echo -e "\t-ppn\t#PE per node (default:2)"
    echo -e "\t-B\tturn on blocking communication (default:disable)"
    echo -e "\t-b\tturn on binding (default:disable)"
    echo -e "\t-l\tnumber of loop (default:100)"
    echo -e "\t-wl\tnumber of warm-up loop (default:10)"
    echo -e "\t-type\ttest type (default:p):\tp -- multi-process, single thread"
    echo -e "\t\t\t\t\ts -- multi-thread (separate context)"
    echo -e "\t\t\t\t\tc -- multi-thread (shared cotnext)"
    echo -e "\t\t\t\t\tk -- multi-thread (shared worker)"
    echo -e "\t-ucx_dev\tUCX_NET_DEVICE (default:mlx5_0:1)"
    echo -e "\t-ucx_tls\tUCX_TLS (default:rc_x)"
    echo -e "\t-v\tprint out commands, do not run commands. (default:disable)"
}

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        -m)
            MSG_RANGE=$VALUE
            ;;
        -w)
            WIN_SIZE=$VALUE
            ;;
        -ppn)
            PE=$VALUE
            ;;
        -B)
            BLOCKING="-B"
            ;;
        -b)
            BINDING="-b"
            ;;
        -type)
            TYPE=$VALUE
            ;;
        -l)
            LOOP=$VALUE
            ;;
        -wl)
            WLOOP=$VALUE
            ;;
        -ucx_dev)
            UCX_DEV=$VALUE
            ;;
        -ucx_tls)
            UCX_TRANS=$VALUE
            ;;
        -v)
            VERBOSE=1
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done

LB="$(cut -d':' -f1 <<<$MSG_RANGE)"
UB="$(cut -d':' -f2 <<<$MSG_RANGE)"
STRIDE="$(cut -d':' -f3 <<<$MSG_RANGE)"

if [[ $TYPE == p ]]
then
    RUN_TYPE="-Dthrds"
    PPR=$PE
    THR=1
elif [[ $TYPE == c ]]
then
    RUN_TYPE="-c"
    PPR=1
    THR=$PE
elif [[ $TYPE == k ]]
then
    RUN_TYPE="-k"
    PPR=1
    THR=$PE
else
    PPR=1
    THR=$PE
fi

NP=$(( 2 * PPR ))

COMMAND="for (( i=$LB;i<=$UB;i*=$STRIDE )) do mpirun --map-by ppr:$PPR:node -np $NP --bind-to none --mca pml ucx --mca osc_base_verbose 0 -x UCX_NET_DEVICES=$UCX_DEV -x UCX_TLS=$UCX_TRANS -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH ./ucx_bench -W $WIN_SIZE -s \$i -n $LOOP -w $WLOOP -t $THR $BINDING $BLOCKING $RUN_TYPE ; done"

echo $COMMAND 2>&1 | tee -a results-$TYPE-$UCX_TRANS-$PE.output

if [[ $VERBOSE == 0 ]]
then
    eval $COMMAND 2>&1 | grep '>' 2>&1 | tee -a results-$TYPE-$UCX_TRANS-$PE.output
fi
