.. SIPPI documentation master file, created by
   sphinx-quickstart on Fri Feb  8 10:00:10 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Setting up SIPPI
================

This section contains information about how to use and control SIPPI,
which requires one to

-  Define the `prior model <#chapprior>`__,

   .. math:: \rho(\mathbf{m})

   ,in form of the prior data structure

   This will allow sampling from a variety of geostatistical models,
   using a unified interface, through ´sippi\_prior.m´

-  Define the `forward model <#chapforward>`__,

   .. math:: g(\mathbf{m})

   , in form of the forward data structure, and the sippi\_forward.m
   m-file

   This will allow solving a 'forward' problem, such as computhing a
   geophysical response from some model.

-  Define the `data and noise model <#chapdata>`__,

   .. math:: L(g(\mathbf{m})

   ,in form of the prior data structure

   This allow describing data and associated uncertainty.

Once these three data structures has been defined, one can solve the
associated inverse problem by `sampling the posterior
distribution </chapSampling/README.md>`__,

.. math:: \sigma(\mathbf{m})


.. toctree::
   :maxdepth: 1
   :caption: Contents:
  
   prior/index
   data/index
   forward/index
   
  
