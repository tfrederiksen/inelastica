[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Build Status](https://travis-ci.org/tfrederiksen/inelastica.svg?branch=master)](https://travis-ci.org/tfrederiksen/inelastica)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/013fe70aa6564ea1bec9df0b3831c834)](https://www.codacy.com/app/brandimarte/inelastica?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tfrederiksen/inelastica&amp;utm_campaign=Badge_Grade)

# Inelastica #

[__Inelastica__][docs] is both the name of this whole Python package
as well as that of an included script to compute inelastic transport characteristics.
[__Inelastica__][docs] is based on the [SIESTA/TranSIESTA][siesta] DFT codes.

The project was initiated around 2003-2005 by Thomas Frederiksen and Magnus Paulsson
while they worked in the group of Mads Brandbyge at the Technical University of Denmark.

## Features ##
[__Inelastica__][docs] contains a number of scripts such as

   - `geom2geom`: Geometry conversion between different file formats
   - `Bandstructures`: Computation of electron and phonon band structures
   - `pyTBT`: A Python version of [tbtrans][siesta] for elastic electron transport
   - `EigenChannels`: Eigenchannel analysis and generation of real-space scattering state wave functions
   - `Phonons`: Vibration modes/frequencies and electron-vibration couplings
   - `Inelastica`: Inelastic transport characteristics (IETS spectroscopy, inelastic shot noise, local heating, etc.)
   - `STM`: Calculation of STM images using the Bardeen approximation

## Dependencies ##
Before installation of Inelastica the following packages are required
   - numpy >= 1.8
   - scipy
   - netCDF4

## Installation ##
Manual installation of Inelastica is performed with the command

    python setup.py install --prefix=<prefix>
    # or
    python setup.py install --home=<my-python-home>

One may also wish to set the following environment variables

    export PYTHONPATH=<my-python-home>/lib/python/
    export PATH=$PATH:<my-python-home>/bin/

## Citations ##
If used to produce scientific contributions please include relevant citations to

    @Article{general-methods,
      Title = {Inelastic transport theory from first principles: Methodology and application to nanoscale devices},
      Author = {Frederiksen, T. and Paulsson, M. and Brandbyge, M. and Jauho, A.-P.},
      Journal = {Phys. Rev. B},
      Year = {2007},
      Number = {20},
      Pages = {205413},
      Volume= {75},
      Doi = {10.1103/PhysRevB.75.205413},
    }
 
    @Article{eigenchannels,
      Title = {Transmission eigenchannels from nonequilibrium Green's functions},
      Author = {Paulsson, M. and Brandbyge, M.},
      Journal = {Phys. Rev. B},
      Year = {2007},
      Number = {11},
      Pages = {115117},
      Volume = {76},
      Doi = {10.1103/PhysRevB.76.115117},
    }

    @Article{stm,
      Title = {Scanning tunneling microscopy current from localized basis orbital density functional theory},
      Author = {Gustafsson, A. and Paulsson, M.},
      Journal = {Phys. Rev. B},
      Year = {2016},
      Pages = {115434},
      Volume = {93},
      Doi = {10.1103/PhysRevB.93.115434},
    }

    @Article{noise,
      Title = {Inelastic shot noise characteristics of nanoscale junctions from first principles},
      Author = {Avriller, R. and Frederiksen, T.},
      Journal = {Phys. Rev. B},
      Year = {2012},
      Pages = {155411},
      Volume = {86},
      Doi = {10.1103/PhysRevB.86.155411},
    }

## Documentation ##
Some documentation may be found [here][docs].

## Contributions, issues and bugs ##
Contributions are highly appreciated.

If you find any bugs please form a [bug report/issue][issues]

If you have a fix please consider adding a [pull request][pulls].

## License ##
The [__Inelastica__][docs] license is [LGPL][lgpl], please see the LICENSE file.

<!---
Links to external and internal sites.
-->
[siesta]: https://launchpad.net/siesta
[issues]: https://github.com/tfrederiksen/inelastica/issues
[pulls]: https://github.com/tfrederiksen/inelastica/pulls
[lgpl]: http://www.gnu.org/licenses/lgpl.html
[docs]: https://tfrederiksen.github.io/inelastica
