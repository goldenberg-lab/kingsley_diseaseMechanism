#!/bin/bash -x

#PBS -N kingsley_Rproject
#PBS -l vmem=60g
#PBS -l walltime=23:50:00
#PBS -joe /hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/scripts/

/hpf/tools/centos6/R/3.2.3/bin/Rscript \
/hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/Rproject/analyse_data.r \
/hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/ \
/hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/grp4_grp3/ \
/hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/grp4_grp3/results/ 0 0
