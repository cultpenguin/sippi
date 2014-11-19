%% simple_prior_models

clear all;close all

im=0;

%% 2D Gaussian, FFTMA
im=im+1;
prior{im}.name='Traditional gaussian (FFTMA)'; % [optional] specifies name to prior
prior{im}.type='FFTMA';                % the type of a priori model
prior{im}.x=[0:1:100];                 % specifies the scales of the 1st (X) dimension
prior{im}.y=[10:1:90];                 % specifies the scales of the 2nd (Y) dimension
prior{im}.m0=2000;                     % the a priori mean value  (default m0=0, if not set)
prior{im}.Cm='200 Sph(30)';     % the a priori covariance/semivariogram model

prior=sippi_prior_init(prior); % initialize the prior model
m=sippi_prior(prior);          % generate a realization from the prior model
sippi_plot_prior(prior,m)      % visualize the realization from the prior
print_mul('prior_example_1_2d_gaussian')

sippi_plot_prior(prior)      % visualize the realization from the prior
print_mul('prior_example_1_2d_gaussian_realizations')


%% 2D Gaussian, VISIM
im=im+1;
prior{im}.name='Traditional gaussian (VISIM)'; % [optional] specifies name to prior
prior{im}.type='VISIM';                % the type of a priori model
prior{im}.x=[0:1:100];                 % specifies the scales of the 1st (X) dimension
prior{im}.y=[10:1:90];                 % specifies the scales of the 2nd (Y) dimension
prior{im}.m0=2000;                     % the a priori mean value  (default m0=0, if not set)
prior{im}.Cm='200 Sph(30)';     % the a priori covariance/semivariogram model

prior=sippi_prior_init(prior); % initialize the prior model
m=sippi_prior(prior);          % generate a realization from the prior model

sippi_plot_prior(prior,m,im)      % visualize the realization from the prior
print_mul('prior_example_2_2d_visim')


%% 2D anisotropic Gaussian, FFTMA
im=im+1;
prior{im}.name='Traditional gaussian (FFTMA)'; % [optional] specifies name to prior
prior{im}.type='FFTMA';                % the type of a priori model
prior{im}.x=[0:1:100];                 % specifies the scales of the 1st (X) dimension
prior{im}.y=[10:1:90];                 % specifies the scales of the 2nd (Y) dimension
prior{im}.m0=2000;                     % the a priori mean value  (default m0=0, if not set)
prior{im}.Cm='200 Sph(30,40,.33)';     % the a priori covariance/semivariogram model

prior=sippi_prior_init(prior); % initialize the prior model
m=sippi_prior(prior);          % generate a realization from the prior model
sippi_plot_prior(prior,m,im)      % visualize the realization from the prior
print_mul('prior_example_3_2d_fftma')


%% 2D SISIM
im=im+1;
prior{im}.name='Sequential Indicator Simulation (SISIM)'; % [optional] specifies name to prior
prior{im}.type='SISIM';                % the type of a priori model
prior{im}.x=[0:1:100];                 % specifies the scales of the 1st (X) dimension
prior{im}.y=[10:1:90];                 % specifies the scales of the 2nd (Y) dimension
prior{im}.Va='1 Gau(30,40,.33)';
prior{im}.marginal_prob=[0.4 0.2 .4];

prior=sippi_prior_init(prior); % initialize the prior model
m=sippi_prior(prior);          % generate a realization from the prior model
sippi_plot_prior(prior,m,im)      % visualize the realization from the prior
print_mul('prior_example_4_2d_sisim')


%% 2D SNESIM
im=im+1;
prior{im}.name='Tranining image base prior (SNESIM)'; % [optional] specifies name to prior
prior{im}.type='SNESIM';                % the type of a priori model
prior{im}.x=[0:1:100];                 % specifies the scales of the 1st (X) dimension
prior{im}.y=[10:1:90];                 % specifies the scales of the 2nd (Y) dimension
prior{im}.S=sgems_get_par('snesim_std'); % Get default SGeMS parameter file, with default training image
prior{im}.scaling=0.5;                 % scaling of the training image
prior{im}.rotation=30;                 % rotation of the training image
prior=sippi_prior_init(prior); % initialize the prior model
m=sippi_prior(prior);          % generate a realization from the prior model
sippi_plot_prior(prior,m,im)      % visualize the realization from the prior
print_mul('prior_example_5_2d_snesim')


