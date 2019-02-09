Sampling the posterior
######################

Once the `prior`, `data`, and `forward` data structures have been
defined, the associated a posteriori probability can be sampled using
[the extended rejection sampler](/chapSampling/chapSampling_rejection.md) and [the extended Metropolis sampler](/chapSampling/chapSampling_metropolis.md).

If the inverse problem is linear and Gaussian it can be solved using [linear least squares](/chapSampling/linear-least-squares.md) inversion.

.. toctree::
   :maxdepth: 1
   :caption: Contents:
  
   sippi_rejection
   sippi_metropolis
   linear-least-squares
  