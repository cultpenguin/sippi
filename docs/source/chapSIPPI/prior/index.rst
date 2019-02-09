.. SIPPI documentation master file, created by
   sphinx-quickstart on Fri Feb  8 10:00:10 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``prior``: The a priori model
=============================

A priori information is defined by the ``prior`` Matlab structure. Any
number of different types of a priori models can be defined. For example
a 1D uniform prior can be defined in ``prior{1}``, and 2D Gaussian prior
can be defined in ``prior{2}``.

Once a prior data structure has been defined (see examples below), a
realization from the prior model can be generated using

::

    m=sippi_prior(prior);

The realization from the prior can be visualized using

::

    sippi_plot_prior(prior,m);

A sample (many realizations) from the prior can be visualized using

::

    m=sippi_plot_prior_sample(prior);

All a priori model types in SIPPI allow to generate a new model in the
vicinity of a current model using

::

    [m_new,prior]=sippi_prior(prior,m);

in such a way that the prior model will be sampled if the process is
repeated (see `Sequential Gibbs Sampling <#sec_seq_gibbs>`__).

Types of a priori models
------------------------

Six types of a priori models are available, and can be selected by
setting the ``type`` in the ``prior`` structure using e.q.
``prior{1}.type='gaussian'``.

The `UNIFORM <uniform.md>`__ type prior specifies an uncorrelated ND
uniform model.

The `GAUSSIAN <gaussian.md>`__ type prior specifies a 1D generalized
Gaussian model.

The `FFTMA <fftma.md>`__ type prior specifies a 1D-3D Gaussian type a
priori model based on the FFT Moving Average method, which is very
efficient for unconditional sampling, and for defining a prior Gaussian
model with variable/uncertain mean, variance, ranges, and rotation.

The `CHOLESKY <cholesky.md>`__ type prior specifies a 1D-3D Gaussian
type a priori model based on Cholesky decomposition of the covariance
model.

The `VISIM <visim.md>`__ type prior model specifies 1D-3D Gaussian
models, utilizing both sequential Gaussian simulation (SGSIM) and direct
sequential simulation (DSSIM) that can be conditioned to data of both
point- and volume support and linear average data.

The `PLURIGAUSSIAN <plurigaussian.md>`__ type prior model specifies
1D-3D pluriGaussian. It is a type if truncated Gaussian model that can
be used for efficient simulation of categorical values.

The `VORONOI <voronoi.md>`__ type prior defines a number of Voronois
cells in a 1D to 3D grid.

The `MPS <mps.md>`__ type prior model specifies a 1D-3D
multiple-point-based statistical prior model, based on the `MPS <#>`__
C++ library. Simulation types includes SNESIM (based on a search tree or
list), ENESIM, and GENESIM (generalized ENESIM).

The `SNESIM <snesim.md>`__ type prior model specifies a 1D-3D
multiple-point-based statistical prior model based on the SNESIM code
from `Stanford/SCRF <#>`__.

The `SNESIM\_STD <snesim_std.md>`__ is similar to the 'SNESIM' type
prior, but is based on `SGEMS <#SGEMS>`__.

The following sectionsdocuments the properties of each type of prior
model.

Examples of using different types of prior models or combining prior
models can be found in the `examples section <#sec_ex_prior>`__.

   
   
.. toctree::
   :maxdepth: 1
   :caption: Contents:

   SequentialGibbs.rst
   uniform.rst
   gaussian.rst	
   fftma.rst
   visim.rst
   cholesky.rst
   plurigaussian.rst
   voronoi.rst
   mps.rst
   snesim.rst
   snesim_std.rst   
   