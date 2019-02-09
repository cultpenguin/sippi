VISIM
-----

::

    im=im+1;
    prior{im}.type='VISIM';
    prior{im}.x=[0:1:100];
    prior{im}.y=[0:1:100];
    prior{im}.m0=10;
    prior{im}.Cm='1 Sph(10)';

.. figure:: ../../figures/prior_visim_2d_gaussian.png
   :alt: 

As with the FFTMA prior model the VISIM prior can make use of a target
distribution. However, if a target distribution is set, the use of the
VISIM prior model will utilize direct sequential simulation, which will
ensure both histogram and covariance reproduction.

Using a target distribution together with the VISIM prior model is
similar to that for the FFTMA prior model. Simply the ``type`` has to be
changed from FFTMA to VISIM:

::

    clear all;close all;
    im=1;
    prior{im}.type='VISIM';
    prior{im}.x=[0:1:40];
    prior{im}.y=[0:1:40];
    prior{im}.m0=10;
    prior{im}.Cm='1 Sph(10)';

    % Create target distribution
    N=10000;
    prob_chan=0.5;
    d1=randn(1,ceil(N*(1-prob_chan)))*.5+8.5;
    d2=randn(1,ceil(N*(prob_chan)))*.5+11.5;
    d_target=[d1(:);d2(:)];
    prior{im}.d_target=d_target;

.. figure:: ../../figures/prior_visim_2d_target.png
   :alt: 

