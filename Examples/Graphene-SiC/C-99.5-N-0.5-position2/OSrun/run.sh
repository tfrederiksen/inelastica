#!/bin/csh
#$ -cwd
#$ -V -S /bin/bash
#$ -N OSrun
#$ -o stdout
#$ -e stdout
#$ -pe smp 24


mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta  < RUN_1.fdf > RUN_1.out 
mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta  < RUN_2.fdf > RUN_2.out 
mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta  < RUN_3.fdf > RUN_3.out 
mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta  < RUN_4.fdf > RUN_4.out 
mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta  < RUN_5.fdf > RUN_5.out 
mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta  < RUN_6.fdf > RUN_6.out 

wait
rm -f fort.*
#end
