#!/bin/bash

#PBS -l walltime=30:00:00
#PBS -l nice=19
#PBS -l nodes=1:ppn=1
#PBS -q batch

cd $PBS_O_WORKDIR
./carreras_sex ${PBS_ARRAYID} 12 36 500 0.01 8000 0.02 0.4 1
