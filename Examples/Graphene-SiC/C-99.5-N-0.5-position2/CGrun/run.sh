#!/bin/csh 
#$ -cwd  
#$ -V -S /bin/bash 
#$ -N CGrun-0.995-VCA 
#$ -o stdout 
#$ -e stdout 
#$ -pe smp 24  
mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta   < ./RUN.fdf > ./RUN.out  
