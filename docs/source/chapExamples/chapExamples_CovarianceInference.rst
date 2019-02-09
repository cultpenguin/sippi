Probilistic covariance/semivariogram indeference
================================================

This chapter documents how to use SIPPI to infer properties of a
covariance/semivariogram model from noisy data (both data of point
support and linear average data can be considered).

To perform probabilistic inference of covariance model parameters one
must

1. define the data and associated uncertainty (if any),

2. define a prior model of the covariance model parameters that one wish
   to infer information about, and

3. define the linear forward operator (only applicable if data are not
   of point support).

The methodology is published in `Hansen et al.
(2015) </bibliography.md>`__.

Specification of covariance model parameters
--------------------------------------------

The following covariance model properties can be defined, that allow
defining an isotropic or an-isotropic covariance model in 1D, 2D, or 3D:

::

    type                % covariance model type (1->Sph, 2->Exp, 3->Gau)
    m0                  % the mean
    sill                % the variance
    nugget_fraction     % percentage of the variance assigned to a Nugget
    range_1             % range in primary direction
    range_2             % range in secondary direction
    range_3             % range in tertiary direction
    ang_1               % first angle of rotation 
    ang_2               % second angle of rotation 
    ang_3               % third angle of rotation 

Inference of a full 1D covariance model requires defining
[type,sill,nugget\_fraction,range\_1].

Inference of a full 2D covariance model requires defining
[type,sill,nugget\_fraction,range\_1,range\_2,ang\_1].

Inference of a full 3D covariance model requires defining
[type,sill,nugget\_fraction,range\_1,range\_2,range\_3,ang\_1,ang\_2,ang\_3].

In order to define which of the covariance model parameters to infer
information about, simply define a prior structure for any of these
parameters, as 1D type SIPPI prior model.

For example, to simple infer information about the range in the primary
direction, with a priori distribution of the range as U[0,3] use

::

    forward.Cm='1 Sph(10)';

    im=1;
    prior{im}.type='uniform';
    prior{im}.name='range_1'; % the 'name' field is used to identify the covariance model parameter!
    prior{im}.min=0;
    prior{im}.max=3;

In this case an ``range_1`` refers to the isotropic range in the
covariance model defined in the ``forward.Cm`` field

If, instead

::

    forward.Cm='1 Sph(10,90,.25)';

then ``range_1`` would refer to the range in the direction of maximum
continuity (90 degrees from North). ``range_2`` will in this case be
fixed.

As described above, the covariance model type can be considered as a
unknown parameter, that can be inferred during inversion. This may pose
some problems as discussed in HCM15.

To infer the covariance model type, a prior 1D structure should be
defined as e.g.

::

    im=1;
    prior{im}.type='uniform';
    prior{im}.name='type'; % 
    prior{im}.min=0;
    prior{im}.max=3;

Any value between 0 and 1 defines a spherical type covariance. Any value
between 1 and 2 defines an exponential type covariance. Any value
between 3 and 3 defines a Gaussian type covariance.

Thus no prior should be defined for the 'type' prior that can provide
values below 0, and above 3. In the case above, all three covariance
model types has the sane a priori probability.

A detailed description of how to parameterize the inverse covariance
model parameter problem, can be found in
`sippi\_forward\_covariance\_inference <#sippi_forward_covariance_inference>`__.

Inferring a 2D covariance model from the Jura data set - Example of point support
---------------------------------------------------------------------------------

The Jura data set (see Goovaerts, 1997) contains a number observations
of different properties in 2D. Below is an example of how to infer
properties of a 2D covariance model from this data set.

A Matlab script implementing the steps below can be found here:
`jura\_covariance\_inference.m <#>`__

Load the Jura data
~~~~~~~~~~~~~~~~~~

Firs the Jura data is loaded.

::

    % jura_covariance_inference
    %
    % Example of inferring properties of a Gaussian model from point data
    %

    %% LOAD THE JURA DATA
    clear all;close all
    [d_prediction,d_transect,d_validation,h_prediction,h_transect,h_validation,x,y,pos_est]=jura;
    ix=1;
    iy=2;
    id=6;

    % get the position of the data
    pos_known=[d_prediction(:,[ix iy])];  

    % perform normal score transformation of tha original data
    [d,o_nscore]=nscore(d_prediction(:,id));
    h_tit=h_prediction{id};

Setting up SIPPI for covariance parameter inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First a SIPPI 'prior' data structure is setup do infer covariance model
parameters for a 2D an-isotropic covariance model. That is, the
``range_1``, ``range_2``, ``ang_1``, and ``nugget_fraction`` are defined
using

::

    im=0;
    % A close to uniform distribution of the range, U[0;3].
    im=im+1;
    prior{im}.type='uniform';
    prior{im}.name='range_1';
    prior{im}.min=0;
    prior{im}.max=3;

    im=im+1;
    prior{im}.type='uniform';
    prior{im}.name='range_2';
    prior{im}.min=0;
    prior{im}.max=3;

    im=im+1;
    prior{im}.type='uniform';
    prior{im}.name='ang_1';
    prior{im}.min=0;
    prior{im}.max=90;

    im=im+1;
    prior{im}.type='uniform';
    prior{im}.name='nugget_fraction';
    prior{im}.min=0;
    prior{im}.max=1;

Thus the a priori information consists of uniform distributions of
ranges between 0 and 3, rotation between 0 and 90, and a nugget fraction
between 0 and 1 is.

Then the data structure is set up, using the Jura data selected above,
while assuming a Gaussian measurement uncertainty with a standard
deviation of 0.1 times the standard deviation of the data:

::

    %% DATA
    data{1}.d_obs=d; % observed data
    data{1}.d_std=0.1*std(d);.4; % uncertainty of observed data (in form of standard deviation of the noise)

Finally the forward structure is setup such that
``sippi_forward_covariance_inference`` allow inference of covariance
model parameters.

In the ``forward`` structure the location of the point data needs to be
given in the ``pos_known`` field, and the initial mean and covariance
needs to be set. Also, the name of the forward function used (in this
case
`sippi\_forward\_covariance\_inference <#sippi_forward_covariance_inference>`__)
must be set. Use e.g.:

::

    %% FORWARD
    forward.forward_function='sippi_forward_covariance_inference';
    forward.point_support=1;
    forward.pos_known=pos_known;

    % initial choice of N(m0,Cm), mean and sill are 0, and 1, due
    % due to normal score
    forward.m0=mean(d);
    forward.Cm=sprintf('%3.1f Sph(2)',var(d));

Now, SIPPI is set up for inference of covariance model parameters. Use
for example the Metropolis sampler to sample the a posterior
distribution over the covariance model parameters using:

::

    options.mcmc.nite=100000;
    options.mcmc.i_plot=1000;
    options.mcmc.i_sample=25;
    options.txt=name;
    [options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)

    sippi_plot_prior(options.txt);
    sippi_plot_posterior(options.txt);

Sampling the posterior provides the following 2D marginal distributions
Note how several areas of high density scatter points (i.e. areas with
high posterior probability) can be found.
