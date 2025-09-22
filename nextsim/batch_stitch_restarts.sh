#! /bin/bash
# Break time period into batches of about 30 mins for multi-node parallelisation.
# (There are 64 cores on each node, so each node can do multiple dates at once.)

function usage {
    echo "Usage:"
    echo "launch_all.sh [START_DATE] [END_DATE] [CORES]"
    echo "    START_DATE: first date to process in yyyymmdd format (default=19910101)"
    echo "    END_DATE: last date to process in yyyymmdd format (default=20241231)"
    echo "    CORES: number of dates to do in parallel (default=32)"
    exit 1
}
[[ "$1" == "-h" ]] || [[ "$1" == "--help" ]] || [[ $# -gt 3 ]] && usage
START_DATE=${1-19910101}
END_DATE=${2-20241231}
CORES=${4-32}
sets=5 # approx 10-20 mins

d1=$START_DATE
while [ $d1 -lt $END_DATE ]
do
    n=$((sets * CORES - 1))
    d2=$(date -d "$d1 + ${n}days" "+%Y%m%d")
    [[ $d2 -gt $END_DATE ]] && d2=$END_DATE
    cmd="sbatch --time=0-00:30:00 --job-name=stitch_${d1}-${d2} slurm.sh $d1 $d2 $CORES"
    echo $cmd
    $cmd
    d1=$(date -d "$d2 + 1 day" "+%Y%m%d")
done
