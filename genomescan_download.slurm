#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_mpu

module load Python/3.8.3-intel-2020a
echo $PYTHONPATH

export PYTHONPATH="/data/antwerpen/205/vsc20587/software/python_lib/lib/python3.8/site-packages/:${PYTHONPATH}"
export gsport=/user/antwerpen/205/vsc20587/data/software/python_lib/bin/gsport


## 36 WGS samples from Allison on quiescence
cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/LeishQ_WGS_allison_20240904/

## download all raw data from the project, including the md5checksum file
$gsport -p 106214 -l ## list all files and see whether recursive needed
$gsport -p 106214 -a -t 2


# calculate the md5sum from the downloaded data
md5sum -c checksums.md5

