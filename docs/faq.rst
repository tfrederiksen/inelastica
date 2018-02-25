.. _faq:

FAQ
===

If you have questions, please send us an email. We will try to help out as much as possible providing that you update this documentation with the help you receive.

EigenChannels
-------------

Device region error
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    ERROR! Too much overlap directly from left-top right
    Make Device Region larger!

This can show up when you are running the Eigenchannel script, either it means just what it saying that Hamiltonian contains overlap between the left and right region, which means you need to make your device region larger or maybe redefine ``TS.TBT.PDOSFrom`` and ``TS.TBT.PDOSTo``.

This error can also show up if there is something wrong with the structure. The atoms have to be organized according to left electrode-center region-right electrode and *pyTBT* also assumes that z-coord of left electrode < z-coord of right electrode, so maybe periodic boundary conditions has to be added to get the coordinates right. A recommended structure in z-direction is

.. code-block:: bash

    ABABabab-mol-ababABAB

where *AB* are unrelaxed and *ab* are relaxed layers.

Wave functions not created
~~~~~~~~~~~~~~~~~~~~~~~~~~

EigenChannels finishes without any errors, but has not written any files with the scattering states (WFs), you have probably encountered the following problem: *The energy where the WFs are requested is not among the energy points in the range specified by (TS.TBT.Emin, TS.TBT.Emax, TS.TBT.NPoints)*. If you do not worry about plotting transmission functions, but just need the WFs at the Fermi energy, use these values:

.. code-block:: bash

    TS.TBT.Emin      0.0 eV
    TS.TBT.Emax      1.0 eV
    TS.TBT.Npoints   1

User description of problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. During the calculation of the toturial(AuH2), I faced following problem. When I type 'EigenChannels EC' in the *TSrun* directory, I got the following message:

   .. code-block:: bash

       Reading geometry from "Chain.XV" file
       SiestaIO.ReadXVFile: Reading Chain.XV
       SiestaIO.ReadXVFile: Reading /home/enixchen/AuH2/L10.00/TSrun/Chain.XV
       Traceback (most recent call last):
        File "/usr/local/bin/EigenChannels", line 4, in <module>
          IE.main()
        File "/usr/local/lib/python2.6/dist-packages/Inelastica/EigenChannelspy",
       line 32, in main
          readbasis()
        File "/usr/local/lib/python2.6/dist-packages/Inelastica/EigenChannelspy",
       line 56, in readbasis
          basis = SIO.BuildBasis(fn[0], general.from_atom, general.to_atom)
        File "/usr/local/lib/python2.6/dist-packages/Inelastica/SiestaIO.py", line
       1160, in BuildBasis
          nn += ions[an].numorb
       KeyError: 79
       When I type 'python Phonon.py' in PHrun directory, following error occured
       in the log file
       Traceback (most recent call last):
        File "PHrun.py", line 15, in <module>>  
        PrintSOrbitals=True)
        File "/usr/local/lib/python2.6/dist-packages/Inelastica/Phonons.py", line
        207, in Analyze
          orbitalIndices,nao =
        GetOrbitalIndices(tree[0][2],atomnumberNotSubstituted)
        File "/usr/local/lib/python2.6/dist-packages/Inelastica/Phonons.py", line
        363, in GetOrbitalIndices
          nao = atomnumber2nao[num]
        KeyError: 79

   **solution**

   `SIESTA`_ will generate the *\*.ion.nc* files, which are needed by `Inelastica <docs_>`_ package, during the calculation if your SIESTA is netCDF supported. Please reinstall your `SIESTA`_ with netCDF supported.


Known problems
--------------

Problems in version 1.1
~~~~~~~~~~~~~~~~~~~~~~~

Using more than one basis set for one type of atom is problematic. The reason is that the *pyTBT* and *EigenChannel* scripts identify the basis for each atom using only the atom number, i.e., not the internal numbering used inside Siesta. This should be solved in `release 1.2 <releases_>`_.

Problems in version 1.0 (but solved in 1.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Eigenchannels* and/or *pyTBT* may crash with confusing error messages. Please note that this behavior could be due to one of the following problems:

* *Eigenchannels* and *pyTBT* needs you to define a device subspace using the ``TS.TBT.PDOSFrom`` / ``TS.TBT.PDOSTo`` keywords.
* Crashes asking for profile library. **Solution:** remove any reference to profile in package directory, reinstall.
* "%include" commands cannot have a comment on the same line. (Assumes the comment is part of the filename ...)

