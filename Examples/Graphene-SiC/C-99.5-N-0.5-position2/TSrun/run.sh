#!/bin/csh 
#$ -cwd  
#$ -V -S /bin/bash 
#$ -N TSrun-VCA 
#$ -o stdout 
#$ -e stdout 
#$ -pe x24 24  
mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta   < ./RUN.fdf > ./RUN.out  
