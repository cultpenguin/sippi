``forward``: The forward model:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The specification of the ``prior`` and ``data`` is intended to be
generic, applicable to any inverse problem considered. The forward
problem, on the other hand, is typically specific for each different
inverse problem.

In order to make use of SIPPI to sample the posterior distribution of an
inverse problem, the solution to the forward problem must be embedded in
a Matlab function with the following input and output arguments:

::

    [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id)

``m`` is a realization of the prior model, and ``prior`` and ``data``
are the Matlab structures defining the prior and the noise model (see
`Prior <#chapprior>`__ and `Data <#chapdata>`__)

``id`` is optional, and can be used to compute the forward response of a
subset of the different types of data available (i.e. ``data{1}``,
``data{2}``,... )

The ``forward`` variable is a Matlab structure that can contain any
information needed to solve the forward problem. Thus, the parameters
for the ``forward`` structure is problem dependent. One option,
``forward.forward_function`` is though generic, and point to the m-file
that implements the forward problem.

The output variable ``d`` is a Matlab structure of the same size of
``data``. Thus, if 4 types of data have been specified, then ``d`` must
also be a structures of size 4.

::

    length(data) == length(d);

Further, ``d{i}`` must refer to an array of the same size as
``data{i}.d_obs``.

An example of an implementation of the forward problem related to a
simple line fitting problem is:

::

    function [d,forward,prior,data]=sippi_forward_linefit(m,forward,prior,data);
        d{1}=forward.x*m{2}+m{1};

This implementation requires that the 'x'-locations, for which the
y-values of the straight line is to be computed, is specified through
``forward.x``. Say some y-data has been observed at locations x=[1,5,8],
with the values [2,4,9], and a standard deviation of 1 specifying the
uncertainty, the forward structure must be set as

::

    forward.forward_function='sippi_forward_linefit';
    forward.x=[1,5,8];

while the data structure will be

::

    data{1}.d_obs=[2 4 9]
    data{1}.d_std=1;

This implementation also requires that the prior model consists of two
1D prior types, such that

::

    m=sippi_prior(prior)

returns the intercept in ``m{1}`` and the gradient in ``m{2}``.

An example of computing the forward response using an intercept of 0,
and a gradients of 2 is then

::

    m{1}=0;
    m{2}=2;
    d=sippi_forward(m,forward)

and the corresponding log-likelihood of m, can be computed using

::

    logL=sippi_likelihood(data,d);

[see more details and examples related to polynomial line fitting at
`polynomial line fitting <#sec_ex_linefit>`__].

The `Examples <#chapExamples>`__ section contains more example of
implementation of different forward problems.

Validating ``prior``, ``data``, and ``forward``
===============================================

A simple way to test the validity of ``prior``, ``data``, and
``forward`` is to test if the following sequence can be evaluated
without errors:

::

    % Generate a realization, m, of the prior model
    m=sippi_prior(prior);
    % Compute the forward response
    d=sippi_forward(m,forward,prior,data);
    % Evaluate the log-likelihood of m
    logL=sippi_likelihood(data,d);


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   sippiforwardlinearm.rst  
   sippiforwardtraveltime.rst
   
	
	