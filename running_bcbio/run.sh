#!/bin/bash
#SBATCH -n 1
#SBATCH -J umccrise_test
#SBATCH -p vccc
#SBATCH --mem=8000
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=80:00:00
source /home/vlad/load_bcbio_prod.sh
bcbio_nextgen.py ../config/bcbio.yaml -n 28 -q vccc -s slurm -t ipython -r 'timelimit=80:00:00' --retries 1 --timeout 120 -p umccrise_test