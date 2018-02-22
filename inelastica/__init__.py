"""
Inelastica: python code for interfacing with SIESTA and TranSIESTA

Provides:
1:	File format conversions for geometries, try: geom2geom --help
2:	Phonon calculations (including e-ph coupling)
3:	Transmission calculations, try: pyTBT --help
4:	Eigenchannel analysis, try: Eigenchannels --help
5:	IETS calculations, try: Inelastica --help
6:	Scripts to set up the above type of calculations.
7:      Generate / read Hamiltonian, Overlap, Green's functions etc

"""

__version__ = filter(str.isdigit, "$Revision: 13 $")
