### DOCUMENTATION

http://sippi.sourceforge.net/

http://sippi.sourceforge.net/sippi.pdf

http://sippi.sourceforge.net/htmldoc/

### INSTALL
To make us of SIPPI, unpack SIPPI.zip to a fodler $SIPPI. 
The start Matlab add the appropriate paths using
	>> addpath $SIPPI
	>> sippi_set_path

Please go to [http://sippi.sourceforge.net] for examples and details on how to use sippi
	
	
### DIRECTORIES :

$SIPPI/
  SIPPI core m-files
  
$SIPPI/plotting
  Directory containing m-files for plotting

  $SIPPI/data
  Currently only contains the data from ARRENAES described in the paper

  $SIPPI/toolboxes
  'fast_marching_kron' -> Fast marching toolbox by Dirk-Jan Kroon. 
         http://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching
  'mGstat' -> geostatistical toolbox for Matlab, by Thomas M Hansen and Knud S Cordua
         http://mgstat.sourceforge.net
  'tomography' -> A few m-files realated to the tomographics examples

$SIPPI/examples/ 
  contains a number of example for using/running SIPPI

$SIPPI/examples/manuscript_examples/prior_tests
  contains m-files used to generate samples from the prior pdf of a number 
  of different prior type models

$SIPPI/examples/manuscript_examples/case_line_fit
  sippi_line_fit.m demonstrates fitting a straight lin, CASE 1 in the manuscript
  
$SIPPI/examples/manuscript_examples/tomography
  contains a number of example or sampling the a posteriori pdf of 
  tomographic inverse problems, CASE 2 in the manuscript

  - using data AM13 (2D)
  sippi_AM13_metropolis_gaussian.m : Metropolis sampling using Gaussian prior
  sippi_AM13_metropolis_bimodal.m : Metropolis sampling using Gaussian prior / bimodal distribution
  sippi_AM13_metropolis_uniform.m : Metropolis sampling using Gaussian prior / uniform distribution
  sippi_AM13_rejection_gaussian.m : Rejection sampling using Gaussian prior
  sippi_AM13_least_squares.m : Least squares inversion using Gaussian prior
  
  - using data AM24 (2D)
  sippi_AM24_gaussian.m : as sippi_AM13_metropolis_gaussian.m but for data set AM24

  - using data AM1234 (3D)
  sippi_AM1234_metropolis_gaussian.m : sippi_AM13_metropolis_gaussian.m but for data set AM1234
  
## Releases history


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
