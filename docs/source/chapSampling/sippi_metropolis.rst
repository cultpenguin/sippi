The extended Metropolis sampler: ``sippi_metropolis.m``
=======================================================

The extended Metropolis algorithm is in general a much more efficient
algorithm (compared to the `rejection
sampler <chapSampling_rejection.md%29>`__ for sampling the a posteriori
probability

The extended Metropolis sampler can be run using

::

    options.mcmc.nite=40000;    % number of iterations, default nite=30000
    options.mcmc.i_sample=50;   % save the current model for every 50 iterations, default, i_sample=500
    options.mcmc.i_plot=1000;   % plot progress of the Metropolis sampler for every 100 iterations
                                % default i_plot=50;
    options.txt='case_line_fit'; % descriptive name appended to output folder name, default txt='';

    [options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)

One can choose to accept all steps in the Metropolis sampler, which
willresult in an algorithm sampling the prior model, using

::

    options.mcmc.accept_all=1; % default [0]

One can choose to accept models that lead to an improvement in the
likelihood, which results in an optimization like algorithm using

::

    options.mcmc.accept_only_improvements=1; % default [0]

See `sippi\_metropolis <#sippi_metropolis>`__ for more details.

Controlling the step length
---------------------------

One optionally, as part of running the `extended Metropolis
sampler <#sec_extendedmetropolis>`__, automatically update the
'step'-length of the `sequential Gibbs sampler <#sec_seq_gibbs>`__ in
order to ensure a specific approximate acceptance ratio of the
Metropolis sampler. See CHM12 for details.

| The default parameters for adjusting the step length, as given below,
  are set in the '`prior.seq\_gibbs <#sec_seq_gibbs_step>`__' structure.
| These parameters will be set the first time 'sippi\_prior' is called
  with the 'prior' structure as output.The default parameters.

::

    prior{m}.seq_gibbs.step_min=0;
    prior{m}.seq_gibbs.step_max=1;
    prior{m}.seq_gibbs.i_update_step=50
    prior{m}.seq_gibbs.i_update_step_max=1000
    prior{m}.seq_gibbs.n_update_history=50
    prior{m}.seq_gibbs.P_target=0.3000

By default, adjustment of the step length, in order to achieve an
acceptance ratio of 0.3 ('prior{m}.seq\_gibbs.P\_target'), will be
performed for every 50 ('prior{m}.seq\_gibbs.i\_update\_step')
iterations, using the acceptance ratio observed in the last 50
('prior{m}.seq\_gibbs.i\_update\_history') iterations.

Adjustment of the step length will be performed only in the first 1000
('prior{m}.seq\_gibbs.i\_update\_step\_max') iterations.

In order to disable automatic adjustment of the step length simply set

::

    prior{m}.seq_gibbs.i_update_step_max=0; % disable automatic step length

Pertubation strategy: Controlling how to perturb the prior (when multiple priors exist)
---------------------------------------------------------------------------------------

When more than one prior structure is defined (e.g. ``prior{1}``,
``prior{2}``,...) the perturbation strategy can be controlled by the
``options.mcmc.pert_strategy`` field.

``options.mcmc.pert_strategy.perturb_all = 0;``: [default] a random
chosen prior is chosen for perturbation. The other prior types are
ignored.

``options.mcmc.pert_strategy.perturb_all = 1;`` Perturb all prior model
at each iteration.

``options.mcmc.pert_strategy.perturb_all = 2;`` Perturb a random
selection of all prior at each iteration.

Further, when ``options.mcmc.pert_strategy.perturb_all = 0;`` the
freqyency with which a prior is perturbed can be set using
``options.mcmc.pert_strategy.i_pert`` and
``options.mcmc.pert_strategy.i_pert_freq``.

For example to perturb only ``prior{1}`` and ``prior{3}``, while
perturbing ``prior{1}`` 4 times more freqeuntly than\ ``prior{3}`` use

::

    options.mcmc.pert_strategy.i_pert = [1,3]; % only perturb prior 1 and 3
    options.mcmc.pert_strategy.i_pert_freq = [2 8]; % perturb prior 3 80% of
                                                    % the time and prior 1 20% % of the time

By default the probability of perturbing a specific prior is uniform. In
case of 3 prior types this is set using:

::

    options.mcmc.pert_strategy.i_pert = [1,2,3]; % only perturb prior 1 and 3
    options.mcmc.pert_strategy.i_pert_freq = [1,1,1]; % same probability of perturbing all prior types


Multiple pertubation strategies
-------------------------------


Several different pertubation strategies can be provided by specifying ``i_pert`` and ``i_pert_freq_freq`` as a cell array. For example to perturb either prior 1 or any of priors 2-9, use 

::

    options.mcmc.pert_strategy.i_pert{1} = [1];
    options.mcmc.pert_strategy.i_pert_freq{1} = [1];
    options.mcmc.pert_strategy.i_pert{2} = [2:9];
    options.mcmc.pert_strategy.i_pert_freq{2} = ones(1,8)./8;
    
By default the probability of choosing a specific strategy is chosen from uniform distribution. A specific frequency of the chosen strategy can be given by.

::

    options.mcmc.pert_strategy.strategy_freq = [.5 .5]; % default is uniform

To choose strategy 2 in 90% of iterations use

::    

    options.mcmc.pert_strategy.strategy_freq = [.1 .9]; % 


The independent extended Metropolis sampler
-------------------------------------------

The 'independent' extended Metropolis sampler, in which each proposed
model is independant of the previously visited model, can be chosen by
forcing the 'step'-length to be 1 (i.e. leading to independant samples
from the prior), using e.g.

::

    % force independent prior sampling
    for ip=1:length(prior);
        prior{ip}.seq_gibbs.step=1;
        prior{ip}.seq_gibbs.i_update_step_max=0;
    end
    % run 'independent' extended Metropolis sampling
    [options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)

Annealing schedule
------------------

Simulated annealing like behavior can be controlled in the
``options.mcmc.anneal`` structure. By default annealing is disabled.

Annealing consist of setting the temperature (similar to scaling the
noise). A temperature foes not affect the exploration. For temperatures
larger than 1, the acceptance ratio increases (the exploration of the
Metropolis sampler increases). For temperatures below 1, the acceptance
ratio decreases (and hence the exploration of the Metropolis sampler).

The temperature is set to ``options.mcmc.anneal.T_begin`` at any
iteration before ``options.mcmc.anneal.i_begin``. The temperature is set
to ``options.mcmc.anneal.T_end`` at any iteration after
``options.mcmc.anneal.i_end``.

In between iteration number ``options.mcmc.anneal.i_start`` and
``options.mcmc.anneal.i_end`` the temperature changes following either
an exponential decay (``options.mcmc.anneal.type='exp'``), or simple
linear interpolation (``options.mcmc.anneal.type='linear'``).

An annealing schedule can be used allow a Metropolis sampler that allow
exploration of more of the model space in the beginning of the chain.
Recall though that the posterior is not sampled until (at least) the
annealing has been ended at iteration, ``options.mcmc.anneal.i_end``, if
the ``options.mcmc.anneal.T_end=1``. This can potentially help not to
get trapped in a local minimum.

To use this type of annealing, where the annealing stops after 10000
iterations, after which the algorithm performs like a regular Metropolis
sampler, use for example

::

    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=10000; %  iteration number when annealing stops

which is equivalent to

::

    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=10000; %  iteration number when annealing stops
    options.mcmc.anneal.T_begin=5; % start temperature
    options.mcmc.anneal.T_end=1; % end temperature

Parallel tempering
------------------

Parallel tempering is implemented according to S13. It is an extension
of the Metropolis algorithm, that start a number of parallel chains of
Metropolis sampling algorithms. Each chain is run with a different
temperature, and the state of each chain is allowed jump between chains
according to some rules that ensure the correct probability density is
sampled. This allow the sampling algorithm to better handle a posterior
distribution with multiple, disconnected, areas of high probability.

The following three setting enable parallel tempering.

::

    % TEMPERING
    options.mcmc.n_chains=3; % set number of chains (def=1, no multiple chains)
    options.mcmc.T=[1 2 3];      % set temperature of chains [1:n_chains]
    options.mcmc.chain_frequency_jump=0.1; % probability allowing a jump between two chains

``options.mcmc.n_chains`` defines the number of chains. If not set only
one chain is used, and the no parallel tempering is performed.

| ``options.mcmc.T`` defines the temperature of each chain. A
  temperature of '1', which is the default, implies no tempering. A
  higher temperature
| allow a chain to be more exploratory.

| ``options.mcmc.chain_frequency_jump`` defines the frequency with which
  a jump from one chain to another is suggested. A value of one means
  that a
| jump is proposed at each iteration, while a value of 0.1 (default)
  means that a jump is only proposed with 10 percentage probability (on
  average
| one in 10 iterations).
