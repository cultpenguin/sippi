% Run all the test for generating a priori realizations


%% FFTMA type 2D Gaussian prior
clear all,close all
prior_reals_fftma

%% FFTMA type 2D Gaussian prior with uncertain covariance model parameters
clear all,close all
prior_reals_fftma_cov

%% FFTMA type 2D Matern type Gaussian prior with uncertain covariance model parameters
clear all,close all
prior_reals_fftma_matern

%% CHOLESKY type 2D Gaussian prior
clear all,close all
prior_reals_cholesky

%% SISIM type (Gaussian indicator simulation) prior model
clear all, close all
prior_reals_sisim

%% SNESIM type prior model
clear all, close all
prior_reals_snesim





%%
clear all, close all
prior_reals_visim


sgems_clean;
visim_clean;