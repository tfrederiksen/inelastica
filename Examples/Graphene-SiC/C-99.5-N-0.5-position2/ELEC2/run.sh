#!/bin/csh 
#$ -cwd  
#$ -V -S /bin/bash 
#$ -N ELEC1-VCA 
#$ -o stdout 
#$ -e stdout 
#$ -pe smp 24  
mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta   < ./RUN.fdf > ./RUN.out  
python postprocessGeneral.py -s ELEC2 -o RUN.out -n 1  
for infile in *.eps; do 
convert $infile ${infile%.*}.png ; 
done 
