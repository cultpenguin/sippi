The rejection sampler: : ``sippi_rejection.m``
==============================================

The rejection sampler provides a simple, and also in many cases
inefficient, approach to sample the posterior distribution.

At each iteration of the rejection sample an independent realization,
m\_pro, of the prior is generated, and the model is accepted as a
realization of the posterior with probability Pacc = L(m\_pro)/L\_max.
It can be initiated using

::

    options.mcmc.nite=400000; % Number of iteration, defaults to 1000
    options.mcmc.i_plot=500; % Number of iteration between visual updates, defaults to 500
    options=sippi_rejection(data,prior,forward,options);

By default the rejection sampler is run assuming a maximum likelihood of
1 (i.e. L\_max = 1). If L\_max is known, then it can be set using in the
``options.Lmax`` or ``options.logLmax`` fields

::

    options.mcmc.Lmax=1e-9;
    options=sippi_rejection(data,prior,forward,options);

or

::

    options.mcmc.logLmax=log(1e-9);
    options=sippi_rejection(data,prior,forward,options);

Alternatively, L\_max can be automatically adjusted to reflect the
maximum likelihood found while running the rejection sampler using

::

    options.mcmc.adaptive_rejection=1
    options=sippi_rejection(data,prior,forward,options);

An alternative to rejection sampling, also utilizing independent
realizations of the prior, that does not require one to set L\_max is
the `independent extended metropolis
sampler <#sec_independentmetropolis>`__, which may be computatinoally
superior to the rejection sampler,
