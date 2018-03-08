'''
Created on 2016/08/09

@author: Emi Minamitani
'''

import os
import subprocess
import re
import argparse
if __name__ == '__main__':

    p=argparse.ArgumentParser(description='PDOS & band plot from siesta output', fromfile_prefix_chars='@')
    p.add_argument("-s", "--sys", metavar='name of system', action='store', type=str, help="system name used in siesta calculation")
    p.add_argument("-n", "--num", nargs='*', metavar='the index of atom to calculate PDOS', action='store', type=int, help="the index of the posirion of the atom to calculate PDOS")
    p.add_argument("-o", "--output", metavar='name of output file', action='store', type=str, help="name of output file from siesta")
    args=p.parse_args()

    print("your system name")
    print(args.sys)

    bandplot='/home/local/siesta-4.0b-485/Util/Bands/new.gnubands'
    pdosplot='/home/local/siesta-4.0b-485/Util/Contrib/APostnikov/fmpdos'
    bandfile=args.sys+'.bands'

    #plotting bandstructure, shift to Fermi level

    banddata=subprocess.check_output([bandplot, '-F', bandfile])

    bandfile=open('bandData.dat', 'w')
    bandfile.write(banddata)
    bandfile.close()

    file=open('bandData.dat', 'r')
    datas=file.readlines()
    for line in datas:
        if (re.search("k_max", line)):
            print(line)
            kinfo=re.split('\s+', line)
            print(kinfo)
            kmin=kinfo[4]
            kmax=kinfo[5]
            print(str(kmin) +","+ str(kmax))

    gnuplotFile=open('bandplot.gp', 'w')
    gnuplotFile.write('set term postscript eps enhanced color "Arial" 25 \n')
    gnuplotFile.write('set output "band.eps" \n')
    gnuplotFile.write('set xrange ['+str(kmin)+' : '+str(kmax)+' ] \n')
    gnuplotFile.write('set yrange [-5:5] \n')
    gnuplotFile.write('plot "bandData.dat" w l noti, \\\n')
    gnuplotFile.write('0 w l lc -1 noti \n')
    gnuplotFile.close()

    subprocess.call(['gnuplot', 'bandplot.gp'])

    #extract the species from PDOS file
    dosFilename=args.sys+'.PDOS'
    dosFile=open(dosFilename, 'r')
    dosline=dosFile.readlines()

    listspecies=[]
    for line in dosline:
        if re.search("species", line):
            species=re.split(r'[="]', line.strip("\n"))
            listspecies.append(species[2])

    list_uniq=list(set(listspecies))

    print(list_uniq)

    #getting the information of Fermi level
    pattern=re.compile(r'^\s*(scf:)\s+(\d+)')
    output=open(args.output, 'r')
    outputdata=output.readlines()
    for line in outputdata:
        if re.match(pattern, line):
            terms=re.split('\s+', line)
            values=[]
            for i in terms:
                if i!="":
                    values.append(i)

    Fermi=float(values[len(values)-2])

    print("Fermi level is "+ str(Fermi))

    #PDOS for respective atomic species
    for atomtype in list_uniq:
        outfile=atomtype+"pdosTot.dat"
        dossettingFile=open("setting.txt", "w")
        dossettingFile.write(dosFilename+"\n")
        dossettingFile.write(outfile+"\n")
        dossettingFile.write(atomtype+"\n")
        dossettingFile.write("0")
        dossettingFile.close()

        os.system(pdosplot+"< setting.txt")

        gnuplotFile=open('pdos.gp', 'w')
        gnuplotFile.write('set term postscript eps enhanced color "Arial" 25 \n')
        gnuplotFile.write('set output  "'+atomtype+'-pdosFromFermilevel.eps" \n')
        gnuplotFile.write('set xrange [-5:5 ] \n')
        gnuplotFile.write('plot  "'+outfile+'" using ($1-'+str(Fermi)+'):2'+' w l noti \n')
        gnuplotFile.close()

        os.system("gnuplot < pdos.gp")

    #PDOS for selected atom

    if (args.num):
        for atom in args.num:
            outfile="atom-"+str(atom)+"-pdos.dat"
            dossettingFile=open("setting.txt", "w")
            dossettingFile.write(dosFilename+"\n")
            dossettingFile.write(outfile+"\n")
            dossettingFile.write(str(atom)+"\n")
            dossettingFile.write("0")
            dossettingFile.close()

            os.system(pdosplot+"< setting.txt")

            gnuplotFile=open('pdos.gp', 'w')
            gnuplotFile.write('set term postscript eps enhanced color "Arial" 25 \n')
            gnuplotFile.write('set output  "'+str(atom)+'-pdos.eps" \n')
            gnuplotFile.write('set xrange [-5:5 ] \n')
            gnuplotFile.write('plot  "'+outfile+'" using ($1-'+str(Fermi)+'):2'+' w l noti \n')
            gnuplotFile.close()

            os.system("gnuplot < pdos.gp")
