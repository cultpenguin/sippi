Cross hole travel time delay computation: The forward problem
=============================================================

A number of different methods for solving the problem of computing the
first arrival travel time of a seismic or electromagnetic wave traveling
between a source in one borehole and a receiver in another borehole has
been implemented in the m-file 'sippi\_forward\_traveltime'.

::

    [d,forward,prior,data]=sippi_forward_traveltime(m,forward,prior,data,id,im)

In order to use this m-file to describe the forward problem specify the
'forward\_function' field in the ``forward`` structure using

::

    forward.forward_function='sippi_forward_traveltime';

In order to use ``sippi_forward_traveltime``, the location of the
sources and receivers must be specified in the ``forward.S`` and
``forward.R``. The number of columns reflect the number of data, and the
number of rows reflect whether data are 2D (2 columns) or 3D (3
columns):

::

    forward.S %  [ndata,ndim]
    forward.R %  [ndata,ndim]

Using for example the data from Arrenæs, the forward geometry can be set
up using

::

    D=load('AM13_data.mat');
    forward.sources=D.S;
    forward.receivers=D.R;

In addition the method used to compute the travel times must also be
given (see below).

In order to use the geometry from the AM13 reference data, and the
Eikonal solution to the wave-equation, the ``forward`` structure can be
defined using

::

    D=load('AM13_data.mat');
    forward.forward_function='sippi_forward_traveltime';
    forward.sources=D.S;
    forward.receivers=D.R;
    forward.type='eikonal';

Ray type forward model (high frequency approximation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ray type models are based on an assumption that the wave propagating
between the source and the receiver has infinitely high frequency.
Therefore the travel time delay is due to the velocity along a ray
connecting the source and receiver.

The linear so-called straight ray approximation, which assumes that the
travel time for a wave traveling between a source and a receiver is due
to the travel time delay along a straight line connecting the source and
receiver, can be chosen using

::

    forward.type='ray';
    forward.linear=1;

The corresponding so-called bended-ray approximation, where the travel
time delay is due to the travel time delay along the fast ray path
connecting a source and a receiver, can be chosen using

::

    forward.type='ray';
    forward.linear=0;

When sippi\_forward\_traveltime has been called once, the associated
forward mapping operator is stored in 'forward.G' such the the forward
problem can simply be solved by calling e.g. 'd{1}=forward.G\*m{1}'

Fat Ray type forward model (finite frequency approximation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fat type model assume that the wave propagating between the source and
the receiver has finite high frequency. This means that the travel time
is sensitive to an area around the raypath, typically defined using the
1st Fresnel zone.

A linear fat ray kernel can be chosen using

::

    forward.type='fat';
    forward.linear=1;
    forward.freq=0.1;

and the corresponding non-linear fat kernel using

::

    bforward.type='fat';
    forward.linear=0;
    forward.freq=0.1;

Note that the center frequency of the propagating wave must also be
provided in the 'forward.freq' field. The smaller the frequency, the
'fatter' the ray kernel.

For 'fat' type forward models we rely on the method described by Jensen,
J. M., Jacobsen, B. H., and Christensen-Dalsgaard, J. (2000).
Sensitivity kernels for time-distance inversion. Solar Physics, 192(1),
231-239

Born type forward model (finite frequency approximation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the Born approximation, considering only first order scattering,
can be chosen using

::

    forward.type='born';
    forward.linear=1;
    forward.freq=0.1;

For a velocity field with small spatial variability one can compute
'born' type kernels (using 'forward.linear=0', but as the spatial
variability increases this is not possible.

For the 'born' type forward model we make use if the method described by
Buursink, M. L., Johnson, T. C., Routh, P. S., and Knoll, M. D. (2008).
Crosshole radar velocity tomography with finite‐frequency Fresnel volume
sensitivities. Geophysical Journal International, 172(1), 1-17.

The eikonal equation (high frequency approximation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The eikonal solution to the wave-equation is a high frequency
approximation, such as the one given above.

However, it is computationally more efficient to solve the eikonal
equation directly, that to used the 'forward.type='ray';' type forward
model.

To choose the eikonal solver to compute travel times use

::

    forward.type='eikonal';

The Accurate Fast Marching Matlab toolbox :
http://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching
is used to solve the Eikonal equation.
