'''
Created on 2016/08/08

@author: Emi Minamitani
'''
import argparse
import os
import re
import shutil

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='electronic properties with VCA setup', fromfile_prefix_chars='@')
    parser.add_argument('fraction', metavar='occupation at C atom, setting of dope', action='store', type=float)
    args = parser.parse_args()

    print(args)
    '''
    optimize the lattice parameter for various VCA set up
    & obtain the band structure, DOS, phonon properties
    '''

    fraction = float(args.fraction)
    os.makedirs('./'+str(fraction))
    os.chdir('./'+str(fraction))

    CpotentialFile = "C"
    NpotentialFile = "N"
    CpotentialFileName = CpotentialFile+".psf"
    NpotentialFileName = NpotentialFile+".psf"
    shutil.copyfile('../'+CpotentialFileName, CpotentialFileName)
    shutil.copyfile('../'+NpotentialFileName, NpotentialFileName)
    shutil.copyfile('../Si.gga.psf', 'Si.gga.psf')
    shutil.copyfile('../postprocess.py', 'postprocess.py')

    #need to modify here
    VCApath = "/home/emi/siesta/siesta-4.1-b2/Util/VCA/mixps"
        ######

    VCAcommand = VCApath+" "+CpotentialFile+"  "+NpotentialFile+" "+str(fraction)
    os.system(VCAcommand)

    formatted = "{0:.5f}".format(fraction)
    #print(formatted)

    synthFile = CpotentialFile+NpotentialFile+"-"+str(formatted)+".synth"
    psfFile = CpotentialFile+NpotentialFile+"-"+str(formatted)+".psf"
    atomTypeinfo = CpotentialFile+NpotentialFile+"-"+str(formatted)

    print("synth filename:"+synthFile)
    print("psf filename:" + psfFile)

    templete = open('../templete.fdf')
    datas = templete.readlines()

    fdfFile = open('4HSiC.fdf', 'w')
    pattern = re.compile(r'^\s*(\d*)\s+(\d*)\s+(C.mpn)')
    for line in datas:

        #line of the C.mpn Chemical Species Label
        match = re.match(pattern, line)
        if match:
            print(match.group(0))
            fdfFile.write("  "+match.group(1)+"  201  "+ atomTypeinfo+"\n")

        else:
            changed = line.replace("C.mpn", atomTypeinfo)
            fdfFile.write(changed)

    #append the data of chemical species and basis

    synth = open(synthFile)
    synthdata = synth.readlines()

    # change the atom index in synthdata to 2
    for i, sdata in enumerate(synthdata):
        if i == 1:
            fdfFile.write("2\n")
        else:
            fdfFile.write(sdata)

    fdfFile.close()

    #also need to modify here
    runFile = open('run.sh', 'w')
    runFile.write("#!/bin/csh \n")
    runFile.write("#$ -cwd \n")
    runFile.write("#$ -V -S /bin/bash \n")
    runFile.write("#$ -N VCA-SiC"+str(fraction)+"\n")
    runFile.write("#$ -o stdout \n")
    runFile.write("#$ -e stdout \n")
    runFile.write("#$ -pe smp 24 \n")
    runFile.write("mpirun  /home/emi/siesta/siesta-4.1-b2/Obj_Transiesta/transiesta   < ./4HSiC.fdf > ./4HSiC.out  \n")
    runFile.write('python postprocessGeneral-newTransiesta.py -s 4HSiC -o 4HSiC.out -n 1 5 \n')
    runFile.write('for infile in *.eps; do \n')
    runFile.write('convert $infile ${infile%.*}.png ; \n')
    runFile.write('done \n')
    runFile.close()

    os.system("qsub run.sh")
