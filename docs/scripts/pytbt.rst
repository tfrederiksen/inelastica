.. _pytbt:

pyTBT
=====

pyTBT is a Python version of TBtrans originally developed by Mads Brandbyge.

usage:
  pyTBT [-h] [-f FN] [-F DEVICEFIRST] [-L DEVICELAST] [-N NPOINTS] [--Emin EMIN] [--Emax EMAX] [-x NK1] [-y NK2] [-a GK1] [-b GK2] [-s] [-j] [-e ETA] [-l ETALEAD] [-d] [--useSigNC] [--NumChan NUMCHAN] [--bulk] [--nobulk] [--scaleSigL SCALESIGL] [--scaleSigR SCALESIGR] [--SpectralCutoff SPECTRALCUTOFF] DestDir

positional arguments:
  DestDir               Destination directory

optional arguments:
  -h, --help            show this help message and exit
  -f FN, --fdf FN       Input fdf-file for TranSIESTA calculation [./RUN.fdf]
  -F DEVICEFIRST, --DeviceFirst DEVICEFIRST
                        First device atom (SIESTA numbering) [TS.TBT.PDOSFrom]
  -L DEVICELAST, --DeviceLast DEVICELAST
                        Last device atom (SIESTA numbering) [TS.TBT.PDOSTo]
  -N NPOINTS, --NPoints NPOINTS
                        Energy points [TS.TBT.NPoints]
  --Emin EMIN           First energy point [TS.TBT.Emin]
  --Emax EMAX           Last energy point [TS.TBT.Emax]
  -x NK1, --Nk1 NK1     k-points Nk1 along a1 [1]
  -y NK2, --Nk2 NK2     k-points Nk2 along a2 [1]
  -a GK1, --Gk1 GK1     Gaussian quadrature k-point sampling for a1 direction
                        (2*GK1+1 points) [0]
  -b GK2, --Gk2 GK2     Gaussian quadrature k-point sampling for a2 direction
                        (2*GK2+1 points) [0]
  -s, --skipsym         Skip inversion (time-reversal) symmetry (i.e., k=-k)
                        that reduces the number of k-point evaluations
  -j, --singlejunction  k-point sample only electrode self-energies
  -e ETA, --eta ETA     Imaginary part added to all energies (device and
                        leads) [1e-06 eV]
  -l ETALEAD, --etaLead ETALEAD
                        Additional imaginary part added ONLY in the leads
                        (surface GF) [0.0 eV]
  -d, --skipDOS         Skip calculation of PDOS
  --useSigNC            Use SigNCfiles
  --NumChan NUMCHAN     Number of eigenchannels [10]
  --bulk                Use bulk in electrodes. The Hamiltonian from the
                        electrode calculation is inserted into the electrode
                        region in the TranSIESTA cell [TS.UseBulkInElectrodes]
  --nobulk              Use only self-energies in the electrodes. The full
                        Hamiltonian of the TranSIESTA cell is used in
                        combination with self-energies for the electrodes
                        [TS.UseBulkInElectrodes]
  --scaleSigL SCALESIGL
                        Scale factor applied to Sigma_L [default=1.0]
  --scaleSigR SCALESIGR
                        Scale factor applied to Sigma_R [default=1.0]
  --SpectralCutoff SPECTRALCUTOFF
                        Cutoff value for SpectralMatrix functions (for
                        ordinary matrix representation set cutoff<=0.0)
                        [default=0.0]
