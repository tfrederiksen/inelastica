.. _eigenchannels:

EigenChannels
=============

Eigenchannels, see Paulsson et al. PRB 76, 115117 (2007)


usage:
  EigenChannels [-h] [-F DEVICEFIRST] [-L DEVICELAST] [-n NUMCHAN] [-B] [-M MOLSTATES] [-r RES] [-w FORMAT] [-e ENERGY] [--eta ETA] [-l ETALEAD] [-f FN] [-s ISPIN] [-x K1] [-y K2] [-u] [--bulk] [--nobulk] [--scaleSigL SCALESIGL] [--scaleSigR SCALESIGR] [--SpectralCutoff SPECTRALCUTOFF] DestDir

positional arguments:
  DestDir               Destination directory

optional arguments:
  -h, --help            show this help message and exit
  -F DEVICEFIRST, --DeviceFirst DEVICEFIRST
                        First device atom (SIESTA numbering) [default:
                        TS.TBT.PDOSFrom]
  -L DEVICELAST, --DeviceLast DEVICELAST
                        Last device atom (SIESTA numbering) [default:
                        TS.TBT.PDOSTo]
  -n NUMCHAN, --NumChan NUMCHAN
                        Number of eigenchannels [default: 4]
  -B, --BothSides       Calculate eigenchannels from both sides [default:
                        False]
  -M MOLSTATES, --MPSH MOLSTATES
                        Calculate eigenstates of the device region Hamiltonian
                        (Molecular Projected Selfconsistent Hamiltonian, MPSH)
                        within [default: +/- 0.0] eV from Ef
  -r RES, --Res RES     Resolution [default: 0.4 Ang]
  -w FORMAT, --format FORMAT
                        Wavefunction format (macu, cube, XSF, or nc) [default:
                        XSF]
  -e ENERGY, --Energy ENERGY
                        Energy where eigenchannel scattering states are
                        evaluated [default: 0.0 eV]
  --eta ETA             Imaginary part added to all energies (device and
                        leads) [default: 1e-06 eV]
  -l ETALEAD, --etaLead ETALEAD
                        Additional imaginary part added ONLY in the leads
                        (surface GF) [default: 0.0 eV]
  -f FN, --fdf FN       Input fdf-file for TranSIESTA calculations [default:
                        ./RUN.fdf]
  -s ISPIN, --iSpin ISPIN
                        Spin channel [default: 0]
  -x K1, --k1 K1        k-point along a1 [default: 0.0]
  -y K2, --k2 K2        k-point along a2 [default: 0.0]
  -u, --useSigNC        Use SigNCfiles [default: False]
  --bulk                Use bulk in electrodes. The Hamiltonian from the
                        electrode calculation is inserted into the electrode
                        region in the TranSIESTA cell [default:
                        TS.UseBulkInElectrodes]
  --nobulk              Use only self-energies in the electrodes. The full
                        Hamiltonian of the TranSIESTA cell is used in
                        combination with self-energies for the electrodes
                        [default: TS.UseBulkInElectrodes]
  --scaleSigL SCALESIGL
                        Scale factor applied to Sigma_L [default=1.0]
  --scaleSigR SCALESIGR
                        Scale factor applied to Sigma_R [default=1.0]
  --SpectralCutoff SPECTRALCUTOFF
                        Cutoff value for SpectralMatrix functions (for
                        ordinary matrix representation set cutoff<=0.0)
                        [default=0.0]
