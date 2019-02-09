Install SIPPI
=============

The latest stable version of SIPPI can be downloaded from http://sippi.sourceforge.net as
`SIPPI.zip <https://sourceforge.net/projects/sippi/files/latest/download?source=files>`__.
This version includes 
 `MGSTAT <https://github.com/cultpenguin/mGstat>`__,
 `VISIM <http://www.sciencedirect.com/science/article/pii/S0098300407001045>`__,
 `SNESIM <https://github.com/SCRFpublic/snesim-standalone>`__, and
 `MPS <https://github.com/ergosimulation/mpslib>`__\ LIB.


---

Unpack ZIPPI.zip somewhere, for example to ``c:\Users\tmh\SIPPI``. Then
setup the Matlab path to point to the appropriate SIPPI directories by
typing:

::

    addpath c:\Users\tmh\SIPPI
    sippi_set_path

| For use on Windows, no other setup should be needed.
| For use on Linux (Ubuntu 16.10), no other setup should be needed.

For use on OS X, Xcode with gcc libraries but be available to run some
of the compiled programs. In addition, the ``DYLD_LIBRARY_PATH`` must be
set to point to the shared libraries neeeded by the compiled programs
for OSX. In MATLAB this can be set using:

::

    setenv('DYLD_LIBRARY_PATH', '/usr/local/bin');

On OSX (and Linux versions that are binary incompatible with Ubuntu
16.10) it may be needed to recompile ``MPSlib`` using:

::

    cd SIPPI/toolboxes/mpslib
    ./configure
    make all

Install SIPPI manually from github
----------------------------------

The current development version of SIPPI (less stable and documented) is
hosted on github and can be downloaded (including MGSTAT, SNESIM, VISIM
and MPS) using git:

::

    cd INSTALL_DIR
    git clone --recursive https://github.com/cultpenguin/sippi.git SIPPI

To update this installation (using mGstat and MSPLIB as submodules) use

::

    git pull --recurse-submodules

Then add a path to SIPPI in Matlab using

::

    addpath INSTALL_DIR/SIPPI
    sippi_set_path

To download SIPPI, mGstat and MPSlib seperately use:

::

    cd INSTALL_DIR
    git clone  https://github.com/cultpenguin/sippi.git SIPPI
    git clone  https://github.com/cultpenguin/mgstat.git mGstat
    git clone  https://github.com/ergosimulation/mpslib.git MPSlib

Then add a path to both SIPPI, mGstat and MPS/matlab using

::

    addpath INSTALL_DIR/mGstat
    addpath INSTALL_DIR/SIPPI
    addpath INSTALL_DIR/MPSlib/matlab
    sippi_set_path

Manual compilation
~~~~~~~~~~~~~~~~~~

| SIPPI (optionally) make use of standalone programs from
| `MPS <http://ergosimulation.github.io/mpslib/>`__
  (`github <https://github.com/ergosimulation/mpslib/>`__),
| SNESIM (`github <https://github.com/SCRFpublic/snesim-standalone>`__)
  and
| VISIM (part of `mGstat <http://mgstat.sourceforge.net/>`__ at
  `github <https://github.com/cultpenguin/mGstat>`__). Pre-compiled
  slef-contained binaries are available for windows and Linux, but for
  use on OS-X one may need to manually these programs.

These programs are needed to make use of the
`MPS </chapSIPPI/prior/mps.md>`__, `VISIM </chapSIPPI/prior/visim.md>`__
and `SNESIM </chapSIPPI/prior/visim.md>`__ type `a priori
models <#prior_types>`__ (and can hence be ignored if not used).

Compiling VISIM
^^^^^^^^^^^^^^^

The source code for VISIM is located in ``mGstat/visim/visim_src``. The
compiled program should be placed in ``mGstat/bin`` and called
``visim``.

VISIM use pre-allocation of memory (set before compilation). This is
defined in ``visim.inc``. If this file is changed, VISIM must be
recomiled.

Compiling SNESIM
^^^^^^^^^^^^^^^^

The source code for SNESIM is available from
`github <https://github.com/SCRFpublic/snesim-standalone>`__.

The compiled program should be placed in ``mGstat/bin`` and called
``snesim.exe`` (windows) ``snesim_glnxa64.exe`` (64 bit Linux)
``snesim_maci64.exe`` (64 but OS X).

Compiling MPS
^^^^^^^^^^^^^

The source code for MPS is available from `github <#>`__. To dowload and
compile MPS use for example

::

    cd INSTALL_DIR
    git clone  https://github.com/ergosimulation/mpslib.git MPSLIB
    cd MPS
    make all

Matlabs path should be updated to include ``INSTALL_DIR/MPSLIB/matlab``.

SGeMS (optional)
~~~~~~~~~~~~~~~~

To make use of the SISIM and SNESIM\_STD type prior models, SGeMS needs
to be available.

Currently only SGeMS version 2.1 (`download <#>`__) for Windows is
supported by SIPPI.

Support for SGeMS will be discontinued after relevase v1.5, as the
`MPS <#>`__ library will be used instead.
