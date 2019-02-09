MPS
---

The 'MPS' type prior make use of the `MPS library <#>`__ for
mulitple-point based simulation. For compilation and installation help
see `Install SIPPI <#InstallSippi>`__. MPS implements the SNESIM (using
both a search tree and a list to stor conditional statistics), and the
generalized ENESIM algoritm. The type of multiple-point algorithm is
define in the method field.

To use the MPS type prior at least the type, dimension, as well as a
training image must be provided:

::

    ip=1; 
    prior{ip}.type='mps';
    prior{ip}.x=1:1:80;
    prior{ip}.y=1:1:80;

A trainin imag emust be set in the 'ti' field, as 1D, 2D, or 3D matrix.
If not set, the classical training image from Strebelle is used,
equivalent to:

::

    prior{ip}.ti=channels;

More examples of traning images are located in the 'mGstat/ti' folder.

MPS provides three different simulation aglrithms, which canbe chosen in
the 'method' field as:

::

    prior{ip}.method='mps_snesim_tree';
    prior{ip}.method='mps_snesim_list';
    prior{ip}.method='mps_genesim';

'mps\_snesim\_tree' is the simulation method selected by default if it
is not specified.

options for MPS
~~~~~~~~~~~~~~~

All options for the MPS type simulation algorithm are be available in
the ``prior{ip}.MPS`` data structure.

For example, to set the number of used multiple grids, set the
``MPS.n_multiple_grids`` as

::

    ip=1; 
    prior{ip}.type='mps';
    prior{ip}.method='mps_snesim';
    prior{ip}.x=0:.1:10;
    prior{ip}.y=0:.1:20;
    [m,prior]=sippi_prior(prior);
    i=0;
    for n_mul_grids=[0:1:4];
        prior{ip}.MPS.rseed=1;
        prior{ip}.MPS.n_multiple_grids=n_mul_grids;
        [m,prior]=sippi_prior(prior);
        i=i+1;subplot(1,5,i);
        imagesc(prior{1}.x,prior{1}.y,m{1});axis image
        title(sprintf('NMG = %d',n_mul_grids));
    end

More details on the use of MPS can be found in the SoftwareX manuscript
that describes MPS.
