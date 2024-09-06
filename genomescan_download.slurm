#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_mpu

module load Python/3.8.3-intel-2020a
echo $PYTHONPATH

export PYTHONPATH="/data/antwerpen/205/vsc20587/software/python_lib/lib/python3.8/site-packages/:${PYTHONPATH}"

# /user/antwerpen/205/vsc20587/data/software/python_lib/bin/gsport
export gsport=/user/antwerpen/205/vsc20587/data/software/python_lib/bin/gsport
# add credentials + token

## go to the correct data directory
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/Gabriel_to_be_renamed_RNAseq_20221025/

# download all raw data from the project, including the md5checksum file
# $gsport -p 105214 -l ## list all files and see whether recursive needed
# $gsport -p 105214 -a -t 28

## go to the correct data directory
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/SINGLE_SuSL_20230116

# download all raw data from the project, including the md5checksum file
# $gsport -p 105328 -l ## list all files and see whether recursive needed
# $gsport -p 105328 -a -t 28

## go to the correct data directory
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/SINGLE_sWGA_20230131/

## download all raw data from the project, including the md5checksum file
# $gsport -p 105360 -l ## list all files and see whether recursive needed
# $gsport -p 105360 -a -t 28


## go to the correct data directory
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/SINGLE_sWGA_20230516/

## download all raw data from the project, including the md5checksum file
# $gsport -p 105561 -l ## list all files and see whether recursive needed
# $gsport -p 105561 -a -t 28


# calculate the md5sum from the downloaded data
# md5sum -c checksums.md5



## go to the correct data directory
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/DIORAPHTE_SuSL_aethiopica_20230814/

## download all raw data from the project, including the md5checksum file
# $gsport -p 105717 -l ## list all files and see whether recursive needed
# $gsport -p 105717 -a -t 28


# calculate the md5sum from the downloaded data
# md5sum -c checksums.md5



## download 8 SureSelect samples from the recent outbreak in Ehtiopia
## go to the correct data directory
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/DIORAPHTE_SuSL_aethiopica_20240416/

## download all raw data from the project, including the md5checksum file
# $gsport -p 106070 -l ## list all files and see whether recursive needed
# $gsport -p 106070 -a -t 28


# calculate the md5sum from the downloaded data
# md5sum -c checksums.md5


## download 8 SureSelect samples from the recent outbreak in Kenia
## go to the correct data directory
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/ProjectUnknown_SuSL_kenia_20240528/

## download all raw data from the project, including the md5checksum file
# $gsport -p 106112 -l ## list all files and see whether recursive needed
# $gsport -p 106112 -a -t 28


# calculate the md5sum from the downloaded data
# md5sum -c checksums.md5


## download 18 SureSelect samples from the experiment of Allison
## go to the correct data directory
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/LeishQ_SuSL_allison_20240528/

## download all raw data from the project, including the md5checksum file
# $gsport -p 106147 -l ## list all files and see whether recursive needed
# $gsport -p 106147 -a -t 28


# calculate the md5sum from the downloaded data
# md5sum -c checksums.md5



## 15 SureSelect samples from SpatialCL - Dioraphte project
# cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/DIORAPHTE_SuSL_aethiopica_20240816/

## download all raw data from the project, including the md5checksum file
# $gsport -p 106211 -l ## list all files and see whether recursive needed
# $gsport -p 106211 -a -t 28


# calculate the md5sum from the downloaded data
# md5sum -c checksums.md5


## 4 SureSelect samples from Adugna - test samples of which 2 are repeated from 
## previous run
cd /user/antwerpen/205/vsc20587/aitg_data/jcdujardin/SomaliLeish_SuSL_aethiopica_20240816/

## download all raw data from the project, including the md5checksum file
$gsport -p 106223 -l ## list all files and see whether recursive needed
$gsport -p 106223 -a -t 2


# calculate the md5sum from the downloaded data
md5sum -c checksums.md5

