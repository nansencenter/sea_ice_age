#! /bin/bash
#SBATCH --account=xxxxx
#SBATCH --job-name=stitch
#SBATCH --output=logs/slurm.%j.log
#SBATCH --time=0-20:00:00
##SBATCH --qos=short
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xxxxx

function usage {
    echo slurm.sh START_DATE END_DATE CORES
    exit 1
}
[[ $# -gt 3 ]] && usage
START_DATE=${1-19910101}
END_DATE=${2-20241231}
CORES=${3-36}

source /cluster/home/timill/pynextsim.apptainer.src

indir=/cluster/work/users/timill/sea_ice_age/restarts/
odir=/cluster/work/users/timill/sea_ice_age/stitched_restarts/

apptainer exec --cleanenv /cluster/projects/nn9878k/sim/apptainer_image_files/nextsim_202311.sif \
    ./batch_stitch_restarts.py $indir $odir $START_DATE $END_DATE \
    --cores $CORES --cores-one-pair=1
