Polymode is a package for the modal analysis of optical fibers, particularly microstructured optical fiber (MOF) or photonic crystal fibers (PCF) with silica or polymer construction.

The performance of such fibres can be calculated from a knowledge of the
effective indices and losses of all or at least a sufficient number of the
modes, to give, for example, local density of states, numerical aperture,
bandwidth or mixing behaviour under various perturbations.

The Polymode software is currently the only open source platform for MOF modal analysis.
It allows easy construction of arbitrary waveguide structures using simple geometric shapes as well as imported and resampled images and utilizes the full rotational symmetry of the waveguide structure.

It aims to:

  * Be flexible, easy and intuitive to use
  * Include multiple modal solvers and modal analysis algorithm
  * Provide a flexible and extensible library for modal calculations
  * Provide a number of tools to optimize and analyze MOF structures

# What can Polymode do? #

Polymode contains pre-defined materials with dispersion profiles.

![http://polymode.googlecode.com/svn/wiki/figures/fp_material_plot.png](http://polymode.googlecode.com/svn/wiki/figures/fp_material_plot.png)

Complicated structures can be quickly built.

![http://polymode.googlecode.com/svn/wiki/figures/fp_wg_structures.png](http://polymode.googlecode.com/svn/wiki/figures/fp_wg_structures.png)

Modes in these structures can be found, plotted and used for further calculation.

![http://polymode.googlecode.com/svn/wiki/figures/fp_mode_examples.png](http://polymode.googlecode.com/svn/wiki/figures/fp_mode_examples.png)

Modes can be tracked over wavelength with adaptive step size to give modal birefringence, dispersion and crossing information.

More details can be found in the [QuickIntroduction](QuickIntroduction.md) wiki.

# Installing #
## Installing on Windows ##
To install Polymode on windows you need to first install Python(x,y) from it's website:
http://pythonxy.com/download.php

Then install Polymode using the windows installer in the [Downloads](http://code.google.com/p/polymode/downloads/list) tab.

## Installing on Linux ##
First install Python, ipython, numpy, scipy and matplotlib. On Ubuntu this can be done issuing the command:
```
sudo aptitude install ipython python-matplotlib python-numpy python-scipy python-setuptools
```

Then download the appropriate egg file for your distribution (32 bit or 64 bit) and install with Setuptools:
```
sudo easy_install {name of egg file}
```

For more details see the GettingStarted wiki.