1D Generalized Gaussian
-----------------------

A 1D generalized Gaussian prior model can be specified using the
'gaussian' type prior model

::

    prior{1}.type='gaussian';

A simple 1D Gaussian distribution with mean 10, and standard deviation
2, can be specified using

::

    ip=1;
    prior{ip}.type='gaussian';
    prior{ip}.m0=10;
    prior{ip}.std=2;

The norm of a generalized Gaussian can be set using the 'norm' field. A
generalized 1D Gaussian with mean 10, standard deviation of 2, and a
norm of 70, can be specified using (The norm is equivalent to the beta
factor referenced in
`Wikipedia:Generalized\_normal\_distribution <#>`__)

::

    ip=2;
    prior{ip}.type='gaussian';
    prior{ip}.m0=10;
    prior{ip}.std=2;
    prior{ip}.norm=70;

A 1D distribution with an arbitrary shape can be defined by setting
``d_target``, which must contain a sample of the distribution that one
would like to replicate. For example, to generate a sample from a
non-symmetric bimodal distribution, one can use e.g.

::

    % Create target distribution
    N=10000;
    prob_chan=0.3;
    d1=randn(1,ceil(N*(1-prob_chan)))*.5+8.5;
    d2=randn(1,ceil(N*(prob_chan)))*.5+11.5;
    d_target=[d1(:);d2(:)];

    % set the target distribution
    ip=3;
    prior{ip}.type='gaussian';
    prior{ip}.d_target=d_target;

The following figure shows the 1D histogram of a sample, consisting of
8000 realizations, generated using

::

    sippi_plot_prior_sample(prior,1:ip,8000);

.. figure:: ../../figures/prior_gaussian_1d.png
   :alt: 

