### DOCUMENTATION

http://sippi.sourceforge.net/

http://sippi.sourceforge.net/sippi.pdf

http://sippi.sourceforge.net/htmldoc/

### INSTALL

#### Latest stable release
The simplest approach to start using SIPPI is to download the latest SIPPI package from http://sippi.sourceforge.net/

Simply download [SIPPI.zip](https://sourceforge.net/projects/sippi/files/latest/download?source=files), and unzip SIPPI.zip to a folder such as $SIPPI.

Then start Matlab and add the appropriate paths using

     >> addpath $SIPPI
     >> sippi_set_path

Please go to (http://sippi.sourceforge.net) for examples and details on how to use SIPPI.	

#### Installation from github
The latest version of SIPPI can be downloaded from github. This is the version the developers use. The documentation is usually slightly outdated compared to the github version. The latest copy of SIPPI, including mGstat and MPSLIB, can be downloaded using (most users of SIPPI would need this):

     cd INSTALL_DIR
     git clone --recursive https://github.com/cultpenguin/sippi.git SIPPI

Then add a path to SIPPI

     >> addpath INSTALL_DIR/SIPPI
     >> sippi_set_path

To update SIPPI, including submodules (mGstat and MPSLIB) use

     git pull --recurse-submodules

     
SIPPI, mGstat and MPSLIB can also be downloaded seperately (a develope weould probably prefer this) from github using 

      cd INSTALL_DIR
      git clone --depth 1 https://github.com/cultpenguin/sippi.git SIPPI
      git clone --depth 1 https://github.com/cultpenguin/mgstat.git SIPPI/toolboxes/mGstat
      git clone --depth 1 https://github.com/ergosimulation/mpslib.git SIPPI/toolboxes/MPSLIB
     

Then add a path to both SIPPI, mGstat in Matlab using:

     >> addpath INSTALL_DIR/SIPPI
     >> addpath INSTALL_DIR/mGstat
     >> addpath INSTALL_DIR/MPSLIB/matlab
     >> sippi_set_path

See also https://www.gitbook.com/book/cultpenguin/sippi for more details on manual installation.

### DIRECTORIES :

**$SIPPI/**
  SIPPI core m-files
  
**$SIPPI/plotting**
  Directory containing m-files for plotting

**$SIPPI/data**
  Currently only contains the data from ARRENAES described in the paper

**$SIPPI/toolboxes'**

+'fast_marching_kron' -> Fast marching toolbox by Dirk-Jan Kroon. 
         http://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching

+'mGstat' -> geostatistical toolbox for Matlab, by Thomas M Hansen and Knud S Cordua
         http://github.com/cultpenguin/mgstat

+'MPSLIB' --> Multiple Point Statistics C++ library, https://github.com/ergosimulation/mpslib

+'tomography' -> A few m-files realated to the cross hole tomographic examples

**$SIPPI/examples/** 
  contains a number of example for using/running SIPPI

**$SIPPI/examples/prior_tests**
  contains m-files used to generate samples from the prior pdf of a number 
  of different prior type models

**$SIPPI/examples/case_line_fit**
  sippi_line_fit.m demonstrates fitting a straight lin, CASE 1 in the SIPPI manuscript
  
**$SIPPI/examples/tomography**
  contains a number of examples of sampling the a posteriori pdf for 
  tomographic inverse problems, CASE 2 in the SIPPI manuscript

* using data AM13 (2D)

sippi_AM13_metropolis_gaussian.m: Metropolis sampling using Gaussian prior

sippi_AM13_metropolis_bimodal.m: Metropolis sampling using Gaussian prior / bimodal distribution

sippi_AM13_metropolis_uniform.m: Metropolis sampling using Gaussian prior / uniform distribution

sippi_AM13_rejection_gaussian.m: Rejection sampling using Gaussian prior

sippi_AM13_least_squares.m: Least squares inversion using Gaussian prior
  
* using data AM24 (2D)

sippi_AM24_gaussian.m: as sippi_AM13_metropolis_gaussian.m but for data set AM24

* using data AM1234 (3D)

sippi_AM1234_metropolis_gaussian.m : sippi_AM13_metropolis_gaussian.m but for data set AM1234

**$SIPPI/examples/covariance_inference**

  - jura_covariance_inference:
  Example of probabilistic covariance model parameter inference following
  Hansen et al., 2015 - A general probabilistic approach for inference of Gaussian model parameters from noisy data of point and volume support. 
  doi:10.1007/s11004-014-9567-5 

  
## Releases history

#v1.5 2016-08-18
Notice that V1.3 was packed without mGstat :/
MPSLIB (https://github.com/ergosimulation/mpslib) added for the first time (MPS based prior sampling, SNESIM and ENESIM type algorithms)

#v1.4 2016-02-01
Notice that V1.3 was packed without mGstat :/
So v1.4 is like V1.4 but repackaged properly with mGstat.
+ sippi_prior_plurigaussian added for the first time.

#v1.3 
Added parallel tempering to sippi_metropolis
Many small bugfixes. 
Testet with Matlab R2015b

#v1.2 
Moving from SVN to GIT (on github)

#1.1.1 (01-12-2014) [rev 272]
Fixed bugh in sippi_get_sample that prevented most sippi_plot_posterior_* algorithms to wokr properly

#1.1 (20-11-2014) [rev 265]
Added consistency checks to sippi_prior_init
Seperated visim, sisim prior types into seperate m-files, sippi_prior_visim, sippi_prior_sisim

#1.03 (22-10-2014) [rev 240]
Removed the use of xcorr from signal processing toolbox
Update sippi_get_sample, and most sippi_plot_posterior_* routines

#1.02 (07-10-2014) [rev 235]
Bug fix (sippi_plot_posterior,  prior{ip}.cax need not be set)

#1.01 (09-07-2014)
LUSIM type a priori model
UNIFORM type a priori model
bug fixes

#1.00 (09-07-2014)
Multiple updates, 
* new figures
* Annealing type schedule for monte Carlo Sampling
* more example
* many bugfixes
* quantifying and accounting for modeling errors

#0.96 (12-03-2014)
Many updates for 
sippi_metropolis
sippi_plot_prior
sippi_plot_posterior


#0.95

#0.94
Updates for handling prior types VISIM/SNESIM/FFTMA/SISIM

#0.93 
Updates for plotting prior and posterior statistics (sippi_plot_prior, sippi_plot_posterior)
Updates for fft_ma for better 1D simulation

# 0.92
Fixed plotting of prior and posterior statistic for a scalar, 1D, 2D, and 3D a priori model types

# 0.90
Initial release
