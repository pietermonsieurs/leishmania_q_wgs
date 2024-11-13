#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=10:40:00
#SBATCH -A ap_itg_mpu

/user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/bin/cnv_depth2cnv.py -c ${cov_file}

## loop over cov file
# for cov_file in *.cov; do sbatch --export=cov_file=${cov_file} /user/antwerpen/205/vsc20587/scratch/leishmania_q_wgs/bin/cnv_depth2cnv_wrapper.slurm; done