
### INSTALL

To make us of SIPPI, unpack SIPPI.zip to a fodler $SIPPI. 
The start Matlab add the appropriate paths using
	>> addpath $SIPPI
	>> sippi_set_path

Please go to [http://sippi.sourceforge.net] for examples and details on how to use sippi
	
	
### DIRECTORIES :

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
