# About SIPPI

SIPPI is a [MATLAB](http://mathworks.com/) toolbox \(compatible with [GNU Octave](https://www.gnu.org/software/octave/)\) that been been developed in order solve probabilistically formulated inverse problems \(Tarantola and Valette, 1982; Tarantola, 2005\) where the solution is the a posteriori probability density


$$
\rho(\mathbf{m}) \ = \ k \ \rho(\mathbf{m}) \ L(g(\mathbf{m})),
$$


where $$g(\mathbf{m})$$ refer to the forward model, $$\rho(\mathbf{m})$$ the a priori model, and $$L(g(\mathbf{m}))$$ the likelihood.

SIPPI allow sampling the a posteriori probability density \(Mosegaard and Tarantola, 1995\) in case the forward model is non-linear, and in case using a combination of a number of widely used geostatistical methods to describe a priori information \(Hansen el al., 2012\).

In order to make use of SIPPI one has to

* [Install](chapInstall/README.md) and setup SIPPI.

* Define [the prior model](/chapSIPPI/prior/README.md), $$\rho(\mathbf{m})$$, in form of the `prior` data structure.

* Define [the forward model](/chapSIPPI/chapSIPPI_forward.md), $$g(\mathbf{m})$$, in form of the `forward` data structure, and the `sippi_forward.m` m-file.

* Define the [data and noise model](/chapSIPPI/chapSIPPI_likelihood.md), i.e. the likelihood $$L(g(\mathbf{m}))$$, in form of the `data` structure.

* Choose a method for [sampling the a posteriori probability density](/chapSampling/README.md) \(i.e. the solution to the inverse problem\).

### Implemented methods and algorithms

A number of different [a priori models](/chapSIPPI/prior/README.md) are available: [UNIFORM](/chapSIPPI/prior/uniform.md), [GAUSSIAN](/chapSIPPI/prior/gaussian.md), [FFTMA](/chapSIPPI/prior/fftma.md), [CHOLESKY](/chapSIPPI/prior/cholesky.md), [VISIM](/chapSIPPI/prior/visim.md), [PLURIGAUSSIAN](/chapSIPPI/prior/plurigaussian.md), [VORONOI](/chapSIPPI/prior/voronoi.md), [MPS](/chapSIPPI/prior/mps.md), [SNESIM](/chapSIPPI/prior/snesim.md).

A number of forward solvers is implented: LINEAR \(linear forward operator\) , TRAVELTIME \(ray, fat, eikonal, born\), GPR\_FW \(full waveform modeling\).

Three methods exist that allow sampling the a posterior probability density: [extended Rejection sampling](/chapSampling/chapSampling_rejection.md), [extended Metropolis sampling](/chapSampling/chapSampling_metropolis.md), and [linear least squares](/chapSampling/linear-least-squares.md).

## Getting started

The best way to learn to use SIPPI is by going through some [examples](/chapExamples/README.md):

* [Lineftting example](/chapExamples/chapExamples_linefitting.md): A simple low-dimensional inverse problem.

* [GPR cross hole tomography](/chapExamples/exGPR/README.md): A more complexe inverse problem illustrating most uses of SIPPI.

## Referencing

Two manuscripts exist describing SIPPI. Part I, is a general introduction on how to setup and use SIPPI. Part II, is an example of using SIPPI to solve cross hole GPR inverse problems \(see  [example](chapExamples/exGPR/README.md)\):

> Hansen, T. M., Cordua, K. S., Looms, M. C., & Mosegaard, K. \(2013\). SIPPI: A Matlab toolbox for sampling the solution to inverse problems with complex prior information: Part 1  — Methodology. Computers & Geosciences, 52, 470-480.  
> DOI:[10.1016/j.cageo.2012.09.004](http://dx.doi.org/10.1016/j.cageo.2012.09.004).
>
> Hansen, T. M., Cordua, K. S., Looms, M. C., & Mosegaard, K. \(2013\). SIPPI: A Matlab toolbox for sampling the solution to inverse problems with complex prior information: Part 2 — Application to crosshole GPR tomography. Computers & Geosciences, 52, 481-492.  
> DOI:[10.1016/j.cageo.2012.09.001](http://dx.doi.org/10.1016/j.cageo.2012.09.001).

The key idea that allow using complex a priori models, referred to as 'sequential Gibbs sampling' is described in detail in

> Hansen, T. M., Cordua, K. S., & Mosegaard, K. \(2012\). Inverse problems with non-trivial priors: Efficient solution through sequential Gibbs sampling. Computational Geosciences, 16\(3\), 593-611.  
> DOI:  
> [doi:10.1007/s10596-011-9271-1](http://dx.doi.oef/doi:10.1007/s10596-011-9271-1)

References to other manuscript considered/used in SIPPI is listed in the [Bibliography](/bibliography.md).

# Acknowledgement

SIPPI make use of other open software projects such as :

* mGstat : [http://mgstat.sourceforge.net](http://mgstat.sourceforge.net)
* MPSLib: [https://github.com/ergosimulation/mpslib](https://github.com/ergosimulation/mpslib)
* VISIM : [http://imgp.nbi.ku.dk/visim.php](http://imgp.nbi.ku.dk/visim.php)
* Accurate Fast Marching Matlab toolbox : [http://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching](http://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching)

Codes and theory has been developed by the [Inverse Modeling and Geostatistics Project](http://imgp.nbi.ku.dk)

