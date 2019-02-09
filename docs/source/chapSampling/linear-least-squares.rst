Linear Least Squares inversion: ``sippi_least_squares.m``
=========================================================

If the prior is defined using a pure (no histogram reproduction)
Gaussian type `prior model </chapSIPPI/prior/README.md>`__, a Gaussian
`likelihood </chapSIPPI/chapSIPPI_likelihood.md>`__/noise for the data,
and a linear forward model, then the a posteriori probability density
will also be Gaussian.

In this case the Gaussian a posterior probability density can be
directly estimated using Linear Least Squares inversion (see e.g.
`Tarantola and Valette (1982) </bibliography.md>`__ or `Tarantola
(2005) </bibliography.md>`__), which is available through
``sippi_least_squares.m``, which can be called using

::

    [m_est,Cm_est,m_reals,,options,data,prior,forward]=sippi_least_squares(data,prior,forward,options);

To compute posterior mean and covariance only use e.g.

::

    [m_est,Cm_est]=sippi_least_squares(data,prior,forward);

A number of realizations from the posterior distribution can also be
computed using

::

    [m_est,Cm_est,m_reals,options]=sippi_least_squares(data,prior,forward);

In this case the computed realizations, as well as all computed data,
will be stored in the folder ``options.txt``, similar to when using
`sippi\_metropolis.m </chapSampling/chapSampling_metropolis.md>`__ and
`sippi\_rejection </chapSampling/chapSampling_rejection.md>`__. Some
figures analyzing the posterior distrbibution can then be generated
using e.g. ``sippi_plot_posterior.m``.

``options.lsq`` contains all the operators that is used for the least
squares inversion (``d0``,\ ``Cd``,\ ``m0``,\ ``Cm``,\ ``G``).
