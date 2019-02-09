AM13 Gaussian, Linear least squares tomography
----------------------------------------------

A Matlab script for the following example is available at
`sippi\_AM13\_least\_squares.m <https://github.com/cultpenguin/sippi/blob/master/examples/case_tomography/sippi_AM13_least_squares.m>`__.

`sippi\_least\_squares.m </chapSampling/linear-least-squares.md>`__
allow solving a linear inverse problem with Gaussian prior and noise
model. The tomographic problem can be considered linear in case any of
the linear forward models are chosen, and the prior parameterized in
slowness.

Load the data

::

    clear all;close all
    D=load('AM13_data.mat');
    txt='AM13';

Define a Gaussian noise model using e.g.:

::

    %% THE DATA
    id=1;
    data{id}.d_obs=D.d_obs;
    data{id}.d_std=D.d_std;

Define a Gaussiain type prior model, using(for example) the FFTMA
method, using slowness (inverse velocity):

::

    im=1;
    prior{im}.type='FFTMA';
    prior{im}.name='Slowness (ns/m)';
    prior{im}.m0=7.0035;
    prior{im}.Va='0.7728 Exp(6)';
    prior{im}.x=[-1:dx:6];
    prior{im}.y=[0:dx:13];
    prior{im}.cax=1./[.18 .1];

Finally, define a linear forward model

::

    forward.forward_function='sippi_forward_traveltime';
    forward.type='ray';forward.linear=1;
    % forward.type='fat';forward.linear=1; % alternative forward model
    % forward.type='born';forward.linear=1; % alternative forward model
    forward.sources=D.S;
    forward.receivers=D.R;
    forward.is_slowness=1; % USE SLOWNESS PARAMETERIZATION

The above represents a linear Gaussian inverse problem. This can be
solved using sampling methods, or it can be solved using `linear least
squares inversion </chapSampling/linear-least-squares.md>`__.

::

    [m_est,Cm_est,m_reals,options]=sippi_least_squares(data,prior,forward,options);
