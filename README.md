# Inelastica #
__Inelastica__ is a Python package for the SIESTA/TranSIESTA DFT codes
as well as a script to compute inelastic transport
(and inelastica is the corresponding repository name).

## Features ##
__Inelastica__ contains a number of scripts such as

   - `geom2geom`: Geometry conversion between different file formats
   - `Bandstructures`: Computation of electron and phonon band structures
   - `pyTBT`: A Python version of tbtrans for elastic electron transport
   - `EigenChannels`: Eigenchannel analysis and generation of real-space scattering state wave functions
   - `Phonons`: Vibration modes/frequencies and electron-vibration couplings
   - `Inelastica': Inelastic transport (IETS)

## Installation ##

### Dependencies ###

## Usage ##
If used to produce scientific contributions please include relevant citations to

    @Article{general-methods,
      Title = {Inelastic transport theory from first principles: Methodology and application to nanoscale devices},
      Author = {Frederiksen, Thomas and Paulsson, Magnus and Brandbyge, Mads and Jauho, Antti-Pekka},
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

## Documentation ##
Some documentation may be found at the Inelastica wiki [here][wiki].


## Contributions, issues and bugs ##
Contributions are highly appreciated.

If you find any bugs please form a [bug report/issue][issues]

If you have a fix please consider adding a [pull request][pulls].

## License ##

The Inelastica license is [LGPL][lgpl], please see the LICENSE file.

<!---
Links to external and internal sites.
-->
[issues]: https://github.com/tfrederiksen/inelastica/issues
[pulls]: https://github.com/tfrederiksen/inelastica/pulls
[lgpl]: http://www.gnu.org/licenses/lgpl.html
[wiki]:  http://dipc.ehu.es/frederiksen/inelastica
