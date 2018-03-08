This directory contains several inputs and scripts for generate the potential for virtual crystal approximation in 4H-SiC.

In virtual crystal approximation, we use the "mixed" pseudo potential to reproduce the doped state.
For the detail, we recommend the document in the following link:
http://personales.unican.es/junqueraj/JavierJunquera_files/Metodos/Pseudos/VCA/VCA.pdf

To generate this mixed potential, we use the mixps program bundled in Util/VCA directory in siesta package.
VCA-mixps-CN.py is a python script to generate mixed pseudo potential for n-doped SiC and calculate band and PDOS by siesta. In this example, we mix C and N potential. This script also automatically submit job to calculate the band structure of bulk 4H-SiC. Thus, please modify several paths and command in the script appropriately to your environment.

If you want to create C 99.5% and N 0.5% pseudo potential, using the following command.

python VCA-mixps-CN.py 0.995

This script generate the pseudopotential named as "CN-0.99500.psf", and running the calculation for 4H-SiC bulk. 
"templete.fdf" is the template file for 4H-SiC structure and calculation setups for siesta calculation.

After the calculation finished, you can get the band structure as included in ref-0.995 directory. 
The band structure is plotted along Gamma-M-K-Gamma line, and you can see that the conduction band crosses the Fermi level around the M point.

Note that the fine band structure with small doping amount like this case depends on the k-mesh. We recommend to check the convergence with respect to k-mesh before you try further calculation using Inelastica.

