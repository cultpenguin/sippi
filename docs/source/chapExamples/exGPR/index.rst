Cross hole tomography
=====================

SIPPI includes a `publically available cross hole GPR data from
Arrenaes <ArrenaesData.md>`__ that is free to be used.

SIPPI also includes the implementation of multiple methods for
`computing the travel time delay between a set of sources and
receivers <TraveltimeForward.md>`__. This allows SIPPI to work on for
example cross hole tomographic forward and inverse problems.

This examples is probably the best way to learnthe capabilities of
SIPPI.

This section contains examples for setting up and running a cross hole
tomographic inversion using SIPPI using the `reference data from
Arrenæs <ArrenaesData.md>`__, different types of a priori and `forward
models <TraveltimeForward.md>`__.

`ArrenaesData.md <ArrenaesData.md>`__

The following section contains a few different examples, and many more
are located in [SIPPI/examples/case\_tomography/].


.. toctree::
   :maxdepth: 1
   :caption: Contents:
  
  
   ArrenaesData.rst  
   TraveltimeForward.rst  
   am13-simple.rst  
   am13-bimodal-a-priori-distribution.rst  
   am13-linear-least-squares.rst  
   
   