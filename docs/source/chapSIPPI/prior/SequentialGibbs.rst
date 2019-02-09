Sequential Gibbs sampling / Conditional Re-sampling
###################################################

All the available types of prior models allow perturbing one realization
of a prior into a new realization of the prior, where the degree of
perturbation can be controlled (from a new independent realization to a
very small change).

This means that a random walk, with an arbitrary 'step-length' can be
performed for any of the a priori types available in SIPPI.

For the a priori types 'FFTMA', 'VISIM', 'CHOLESKY', 'SISIM', 'SNESIM',
sequential Gibbs sampling HCM12 is applied. Sequential Gibbs is in
essence a type of conditional re-simulation. From a current realization
of a prior model, a number of model parameters are discarded and treated
as unknown. The unknown model parameters are then re-simulated
conditional to the known model parameters.

In order to generate a new realization 'm2' in the vicinity of the
realization 'm1' use

::

    [m1,prior]=sippi_prior(prior);
    [m2,prior]=sippi_prior(prior,m1);

If this process is iterated, then a random walk in the space of a priori
acceptable models will be perform. Moreover, the collection of
realization obtained in this way will represent a sample from prior
distribution.

Note that in order to use sequential Gibbs sampling ``prior`` must be
given both as an input variable, and as an (possibly update) output
variable.

Controlling sequential Gibbs sampling / Conditional Re-sampling
---------------------------------------------------------------

All properties related to sequential Gibbs sampling can be set in the
'seq\_gibbs' structure (which will be avaiable the first time
`sippi\_prior <#sippi_prior>`__ is called, or if
`sippi\_prior\_init <#sippi_prior_init>`__ is called), for the
individual prior models.

The step-length (i.e. the degree of perturbation) is determined by the
prior{m}.seq\_gibbs.step\` parameter.

For the 'uniform' and 'gaussian' type a priori models a step-length
closer to 0 zeros imples a 'shorter' step, while a step-length close to
1, implies a 'longer' step-length. A step length of 1, will generate a
new independent realization of the prior, while a step length of 0, will
return the same realization of the prior

::

    prior{m}.seq_gibbs.step=.1;
    [m2,prior]=sippi_prior(prior,m1);

For the 'FFTMA', 'VISIM', 'CHOLESKY', 'SISIM', and 'SNESIM' type a
priori models two types (defined in the ``prior{m}.seq_gibbs.type``
variable).

The default 'type' is 2, defined as

::

    prior{m}.seq_gibbs.step=1;
    prior{m}.seq_gibbs.type=2;

where the step length defines the percentage of the of model parameters
(selected at random) defined in ``prior{im}`` is conditionally
re-sampled. Thus, a step-length closer to 0 zeros imples a 'shorter'
step, while a step-length close to 1, implies a 'longer' step-length.

If ``prior{m}.seq_gibbs.step=1``, then ``prior{m}.seq_gibbs.step``
defines the size of a square rectangle/cube which is to be conditionally
re-simulated using sequential Gibbs sampling.
