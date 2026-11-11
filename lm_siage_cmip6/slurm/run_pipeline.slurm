#!/bin/bash
#SBATCH --job-name=cmip6-siage
#SBATCH --account=YOUR_NIRD_ACCOUNT
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/cmip6-siage-%j.out
#SBATCH --error=logs/cmip6-siage-%j.err

set -euo pipefail

CONFIG=${CONFIG:-config/noresm2_mm_historical_r1i1p1f1.yaml}
CONDA_ENV=${CONDA_ENV:-lm-siage-cmip6}

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

export MPLBACKEND=Agg
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mkdir -p logs
python scripts/run_pipeline.py --config "$CONFIG"
