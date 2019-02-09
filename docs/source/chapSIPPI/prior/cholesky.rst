CHOLESKY - 3D Gaussian model
----------------------------

The CHOLESKY type prior utilizes Cholesky decomposition of the
covariance in order to generate realizations of a Gaussian random field.
The CHOLESKY type prior needs a full description of the covariance
model, which will be of size [nxyz\*nxyz\*nxyz], unlike using the
`FFTMA <#prior_fftma>`__ type prior model that only needs a
specification of an isotropic covariance models of size ``[1,nxyz]``.
Hence, the CHOLESKY type prior is much more demanding on memory, and
CPU. However, the CHOLESKY type prior can be used to sample from any
covariance model, also non-stationary covariance model.

The CHOLESKY model is can be defined almost identically to the
`FFTMA <#prior_fftma>`__ type prior model. As an example:

::

    im=1;
    prior{im}.type='CHOLESKY';
    prior{im}.x=[0:2:100];
    prior{im}.y=[0:2:100];
    prior{im}.m0=10;
    prior{im}.Cm='1 Sph(10)';

.. figure:: ../../figures/prior_cholesky_2d.png
   :alt: 

the use of ``d_target`` to specify target distributions is also
possible, using the same style as for the `FFTMA <#prior_fftma>`__ type
prior.

Be warned that the 'cholesky' type prior model is much more memory
demanding than the 'fftma' and 'visim' type prior models, as a full
[nxyz\*nxyz] covariance model needs to setup (and inverted). Thus, the
'cholesky' type prior is mostly applicable when the number of model
parameters (nx\*ny\*nx) is small.
