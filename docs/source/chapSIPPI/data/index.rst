``data``: Data and data uncertainties/noise
===========================================

``data`` is a Matlab structure that defines any number of data and the
associated uncertainty/noise model.

``data{1}`` defines the first data set (which must always be defined),
and any number of additional data sets can be defined in ``data{2}``,
``data{3}``, ...

This allows to consider for example seismic data in ``data{1}``, and
electromagnetic data in ``data{2}``.

For each set of data, a Gaussian noise model (both correlated and
uncorrelated) can be specified. The noise model for different data types
(e.g. ``data{1}`` and ``data{2}`` are independent).

Once the noise model has been defined, the log-likelihood related to any
model, ``m``, with the corresponding `forward
response <#chapforward>`__, ``d``, can be computed using

::

    [d,forward,prior,data]=sippi_forward(m,forward,prior,data)
    logL=sippi_likelihood(data,d)

where ``d`` is the output of `sippi\_forward <#sippi_forward>`__.

The specification of the noise model can be divided into a description
of the `measurement noise <#sec_meas_noise_gauss>`__ (mandatory) and the
`modeling error <#sec_model_noise_gauss>`__ (optional).

Gaussian measurement noise
--------------------------

Uncorrelated Gaussian measurement noise
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To define a set of observed data, [0,1,2], with an associated
uncorrelated uncertainty defined by a Gaussian model with mean 0 and
standard deviation 2, use

::

    data{1}.d_obs=[0 1 2]';
    data{1}.d_std=[2 2 2]';

which is equivalent to (as the noise model for each data is the same,
and independent)

::

    data{1}.d_obs=[0 1 2]';
    data{1}.d_std=2;

One can also choose to define the uncertainty using a variance as
opposed to the standard deviation

::

    data{1}.d_obs=[0 1 2]';
    data{1}.d_var=4;

Correlated Gaussian measurement noise
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Correlated Gaussian measurement uncertainty can be specified using the
``Cd`` field, as for example

::

    data{1}.Cd=[4 1 0 ; 1 4 1 ; 0 1 4];

Note that ``data{1}.Cd`` must be of size [NDxND], where ND is the number
of data in ``data{1}.d_obs``.

Gaussian modeling error
-----------------------

The modeling error refers to errors caused by using for example an
imperfect forward model, see HCM14.

A Gaussian model of the modeling error can be specified by the mean,
``dt``, and the covariance, ``Ct``.

For example

::

    data{1}.dt=[0 0 0];
    data{1}.Ct=[4 4 4; 4 4 4; 4 4 4];

is equivalent to

::

    data{1}.Ct=4

which implies a zero mean modeling error with a covariance model where
all model parameters has a covariance of 4.

`sippi\_compute\_modelization\_forward\_error <#sippi_compute_modelization_forward_error>`__
can be used to estimate the modeling error related to using an
approximate forward model. See the `tomography
example <#sec_ex_tomography>`__, for an `example of accounting for
correlated modeling errors <#AM13_gaussian_modeling_error>`__, following
HCM14.
