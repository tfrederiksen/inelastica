#!/usr/bin/env python

import Inelastica.EigenChannels as EC
import Inelastica.pyTBT as TBT
import Inelastica.Phonons as P
import Inelastica.iets as IETS
import Inelastica.info as info

ver = info.label

# Loop over the three orientations of the carbon chain
for d in ['A1', 'A2', 'A3']:

    # Eigenchannels
    my_argv = '-f %s/TSrun/RUN.fdf %s/%s/ECscript'%(d, d, ver)
    opts = EC.GetOptions(my_argv)
    EC.main(opts)

    # pyTBT
    my_argv = '--Emin=-5.0 --Emax=5.0 -N 11 -f %s/TSrun/RUN.fdf %s/%s/TBTscript'%(d, d, ver)
    opts = TBT.GetOptions(my_argv)
    TBT.main(opts)

    # Phonons
    my_argv = '--FCwildcard=%s/FC* --OSdir=%s/OSrun'%(d, d)
    my_argv += ' --FCfirst=9 --FClast=9 --DeviceFirst=8 --DeviceLast=13 -c %s/%s/PHscript'%(d, ver)
    opts = P.GetOptions(my_argv)
    P.main(opts)

    # IETS
    my_argv = '-F 8 -L 13 -p %s/%s/PHscript/Output.nc -f %s/TSrun/RUN.fdf %s/%s/INscript'%(d, ver, d, d, ver)
    opts = IETS.GetOptions(my_argv)
    IETS.main(opts)
