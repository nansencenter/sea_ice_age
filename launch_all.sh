#! /bin/bash
# Break time period into batches of about 30 mins for multi-node parallelisation.
# (There are 64 cores on each node, so each node can do multiple dates at once.)

cores=36
sets=5 # approx 10-20 mins
start_date=19930101
end_date=20241231

d1=$start_date
while [ $d1 -lt $end_date ]
do
    n=$((sets * cores - 1))
    d2=$(date -d "$d1 + ${n}days" "+%Y%m%d")
    [[ $d2 -gt $end_date ]] && d2=$end_date
    echo sbatch --time="0-00:30:00" --job-name="stitch_${d1}-${d2}" slurm.sh $d1 $d2 $cores
    sbatch --time="0-00:30:00" --job-name="stitch_${d1}-${d2}" slurm.sh $d1 $d2 $cores
    d1=$(date -d "$d2 + 1 day" "+%Y%m%d")
done
