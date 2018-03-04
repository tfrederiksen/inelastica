.. _inelastica:

Inelastica
==========

Inelastica script calculates and writes LOE quantities in ascii (Systemlabel.IN) and NetCDF (Systemlabel.IN.nc)

usage:
  Inelastica [-h] [-n NUMCHAN] [-F DEVICEFIRST] [-L DEVICELAST] [-e ENERGY] [--eta ETA] [-f FN] [-s ISPIN] [-x K1] [-y K2] [-p PHONONNETCDF] [-t TEMP] [-b BIASPOINTS] [-v MAXBIAS] [-c MODECUTOFF] [-V VRMS] [-H] [-d PHEXTDAMP] [-u] [-l ETALEAD] [--SpectralCutoff SPECTRALCUTOFF] [--bulk] [--nobulk] [--scaleSigL SCALESIGL] [--scaleSigR SCALESIGR] [--LOEscale LOESCALE] [--VfracL VFRACL]  DestDir

positional arguments:
  DestDir               Destination directory

optional arguments:
  -h, --help            show this help message and exit
  -n NUMCHAN, --NumChan NUMCHAN
                        Number of eigenchannels [default: 4]
  -F DEVICEFIRST, --DeviceFirst DEVICEFIRST
                        First device atom (SIESTA numbering) [TS.TBT.PDOSFrom]
  -L DEVICELAST, --DeviceLast DEVICELAST
                        Last device atom (SIESTA numbering) [TS.TBT.PDOSTo]
  -e ENERGY, --Energy ENERGY
                        Energy reference where Greens functions etc are
                        evaluated [default: 0.0 eV]
  --eta ETA             Tiny imag. part in Greens functions etc. [default:
                        1e-06 eV]
  -f FN, --fdf FN       Input fdf-file for TranSIESTA calculations [default:
                        ./RUN.fdf]
  -s ISPIN, --iSpin ISPIN
                        Spin channel [default: 0]
  -x K1, --k1 K1        k-point along a1 [default: 0.0]
  -y K2, --k2 K2        k-point along a2 [default: 0.0]
  -p PHONONNETCDF, --PhononNetCDF PHONONNETCDF
                        Electron-phonon coupling NetCDF [default: Output.nc]
  -t TEMP, --Temp TEMP  Temperature [default: 4.2 K]
  -b BIASPOINTS, --BiasPoints BIASPOINTS
                        Number of bias points [default: 801]
  -v MAXBIAS, --MaxBias MAXBIAS
                        Sets the IETS bias range (-MaxBias to MaxBias)
                        [default: 0.4 V]
  -c MODECUTOFF, --ModeCutoff MODECUTOFF
                        Ignore phonon modes with lower hw [default: 0.0025 eV]
  -V VRMS, --Vrms VRMS  Lock in amplifier broadening [default: 0.005 V]
  -H, --Heating         Include heating of vibrational modes [default: False]
  -d PHEXTDAMP, --PhExtDamp PHEXTDAMP
                        External damping [default: 1e-15 (?) TODO check unit!]
  -u, --useSigNC        Use SigNCfiles [default: False]
  -l ETALEAD, --etaLead ETALEAD
                        Additional imaginary part added ONLY in the leads
                        (surface GF) [default: 0.0 eV]
  --SpectralCutoff SPECTRALCUTOFF
                        Cutoff value for SpectralMatrix functions (for
                        ordinary matrix representation set cutoff<=0.0)
                        [default: 1e-08]
  --bulk                Use bulk in electrodes. The Hamiltonian from the
                        electrode calculation is inserted into the electrode
                        region in the TranSIESTA cell [TS.UseBulkInElectrodes]
  --nobulk              Use only self-energies in the electrodes. The full
                        Hamiltonian of the TranSIESTA cell is used in
                        combination with self-energies for the electrodes
                        [TS.UseBulkInElectrodes]
  --scaleSigL SCALESIGL
                        Scale factor applied to Sigma_L [default: 1.0]
  --scaleSigR SCALESIGR
                        Scale factor applied to Sigma_R [default: 1.0]
  --LOEscale LOESCALE   Scale factor to interpolate between LOE-WBA (0.0) and
                        generalized LOE (1.0), see PRB 89, 081405(R) (2014)
                        [default: 1.0]
  --VfracL VFRACL       Voltage fraction over the left-center interface
                        [default: 0.5]
