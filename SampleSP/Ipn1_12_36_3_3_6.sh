#!/bin/bash

#PBS -l walltime=248:00:00
#PBS -l nice=19
#PBS -l nodes=1:ppn=1
#PBS -q linux-spool

cd $PBS_O_WORKDIR
./rwn1sp ${PBS_ARRAYID} 
