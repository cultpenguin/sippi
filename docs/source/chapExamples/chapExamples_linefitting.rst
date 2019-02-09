Polynomial line fitting
=======================

An simple of application of SIPPI is to perform linefitting, as a
probabilistic inverse problem.

Here follows simple polynomial (of order 0, 1 or 2) line-fitting is
considered. Example m-files and data can be found in the
``SIPPI/examples/case_linefit`` folder.

First, the forward problem is defined. Then examples of stochastic
inversion using SIPPI is demonstrated using a a synthetic data set.

The challenge
=============

Assume that some observed data ``d_obs`` are available, as function of a
corresponding set of model parameters ``x``. Assume alse that the
observed data are contaminated by Gaussian noise, with mean 0, and
standard devaiation 10.

::

::

    cd SIPPI/examples/case_linefit
    load sippi_linefit_data
    errorbar(x,d_obs,d_std)

or

::

::

    x=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20];
    d_obs = [ -36.5  -10.8  -24.0  -22.5  -16.1  -10.7   -9.2   -1.6   -6.4   -0.7 -18.5];
    d_std = ones(size(d_obs)).*10;
    errorbar(x,d_obs,d_std)

.. figure:: ../figures/sippi\_linefit\_data\_11.png :alt:

The problem is now to infer information about the 3 parameters defining
a 2nd order polynomium, given the data above.

Solving the problem using SIPPI
===============================

In roder to use SIPPI, the ``data``, ``prior``, and ``forward`` data
strutures must be setup, and the a way of solving the forward problem
must be defined.

The data
--------

The observed data, as well as the associated uncertainty id defined in
the ``sippi`` structure

::

::

    data{1}.d_obs=d_obs;
    data{1}.d_std=d_std;

The prior
---------

A prior distributions of tree parameters, reflecting coefficient of the
polynomium, can be defined as for example;

::

::

    %% setup the prior model
    % the intercept
    im=1;
    prior{im}.type='gaussian';
    prior{im}.name='intercept';
    prior{im}.m0=0;
    prior{im}.std=30;

    % 1st order, the gradient
    im=2;
    prior{im}.type='gaussian';
    prior{im}.name='gradient';
    prior{im}.m0=0;
    prior{im}.std=4;
    prior{im}.norm=80;

    % 2nd order
    im=3;
    prior{im}.type='uniform';
    prior{im}.name='2nd';
    prior{im}.min=-.3;
    prior{im}.max=.3;

Note that each model parameter is associated with it ownm independent,
distritution. Here two Gaussian (N(0,30^2), and N(0,4^2)) as well as a
uniform (U[-0.3 0.3]) is assumed.

A sample of the prior can be generated using

::

::

    m=sippi_prior(prior)

    m =

      1×3 cell array

        [-19.4704]    [3.1968]    [-0.1655]

The forward problem
-------------------

Finally, a way to solvge the forward problem must be implemented, that
takes as input, at least ``m`` (a realization of the pror) and
``forward`` (a matlab structure that contains any informaion neede to
solve the forward problem), and that produces an output ``d``, where
``d{1}`` is of exactly the same size and type as ´data{1}.d\_obs´.

One way to implement this as the m-file
``[sippi_forward_linefit.m](#sippi_forward_linefit)``:

::

::

    % sippi_forward_linefit Line fit forward solver for SIPPI 
    %
    % [d,forward,prior,data]=sippi_forward_linefit(m,forward,prior,data);
    %
    function [d,forward,prior,data]=sippi_forward_linefit(m,forward,prior,data);

    if length(m)==1;
        d{1}=forward.x.*0 + m{1};
    elseif length(m)==2;
        d{1}=forward.x*m{2}+m{1};
    else
        d{1}=forward.x.^2*m{3}+forward.x*m{2}+m{1};
    end

Here ``forward.x`` must be an array of the x-locations, for which the
d-values will be computed.

Note that the prior must be defined such that ``prior{1}`` refer to the
intercept, ``prior{2}`` to the gradient, and ``prior{3}`` to the 2nd
order polynomial coefficient.

If only one prior type is defined then the forward response will just be
a constant, and if two prior types are defined, then the forward
response will be a straight line.

Having implemented the m-file that solves the forward problem in the
style reqiured by SIPPI, the forward can be setup using

::

::

    %% Setup the forward model in the 'forward' structure
    forward.x=x
    forward.forward_function='sippi_forward_linefit';

Evaluate the prior, data, and forward
-------------------------------------

A simple way find problems related to how ``prior``, ``data``,
``forward``, and ``sippi_forward_linefit`` has been setup, correctly is
to test wether he follwing three lines can be executed without errors:

::

::

    m=sippi_prior(prior);
    d=sippi_forward(m,forward);
    logL=sipppi_likelihood(d,data);

Sampling the a posterior distribution
=====================================

Information about the model parameters can be inferred by running the
``extended Metropolis sampler <#metropolis>``\ \_\_ or the
``rejection sampler <#rejection>``\ \_\_ using

Using the Metropolis sampler
----------------------------

The ``extended Metropolis sampler <chapSampling_metropolis.md>``\ \_\_
can be setup and run using:

::

::

    options.mcmc.nite=40000;  % Run for 40000 iterations
    options.mcmc.i_sample=100; % Save every 100th visited model to disc
    options.mcmc.i_plot=5000; % Plot the progress information for every 2500 iterations
    options.txt='case_line_fit_2nd_order'; % descriptive name for the output folder

    [options]=sippi_metropolis(data,prior,forward,options);

Generic statistics about the posterior can be plotted using.

::

::

    % plot posterior statistics, such as 1D and 2D marginals from the prior and posterior distributions
    sippi_plot_prior_sample(options.txt);
    sippi_plot_posterior(options.txt);
    20140521_1644_sippi_metropolis_case_line_fit_2nd_order_m1_3_posterior_sample.png

The figure below show the prior and posterior distribution of the 3
model parameters, as well as the reference values (used to generate the
synthetic data set, in green)

.. figure:: /assets/sippi\_linefit\_data\_11\_postmarg.png :alt:

The figure below plots forward response related to the obtained sample
of the posterior distribution over the model parameters (gray), as well
as the observed data (black), and the noise free reference data obtained
from the reference set of model parameters. \|image0\|

Using the rejection sampler
---------------------------

In a similar manner the
``rejection sampler <chapSampling_rejection.md>``\ \_\_ can be setup and
run using

::

::

    options.mcmc.adaptive_rejection=1; % automatically adjust the normalizing likelihood
    options.mcmc.nite=100000;
    options=sippi_rejection(data,prior,forward,options);

.. \|image0\| image:: /assets/sippi\_linefit\_data\_11\_post.png
