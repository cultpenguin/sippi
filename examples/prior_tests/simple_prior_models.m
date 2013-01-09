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
sippi_plot_model(prior,m)      % visualize the realization from the prior
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

sippi_plot_model(prior,m,im)      % visualize the realization from the prior
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
sippi_plot_model(prior,m,im)      % visualize the realization from the prior
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
sippi_plot_model(prior,m,im)      % visualize the realization from the prior
print_mul('prior_example_4_2d_sisim')


%% 2D SNESIM
im=im+1;
prior{im}.name='Sequential Indicatot Simulation (SISIM)'; % [optional] specifies name to prior
prior{im}.type='SNESIM';                % the type of a priori model
prior{im}.x=[0:1:100];                 % specifies the scales of the 1st (X) dimension
prior{im}.y=[10:1:90];                 % specifies the scales of the 2nd (Y) dimension
prior{im}.S=sgems_get_par('snesim_std'); % Get default SGeMS parameter file, with default training image
prior{im}.scaling=0.5;                 % scaling of the training image
prior{im}.rotation=30;                 % rotation of the training image
prior=sippi_prior_init(prior); % initialize the prior model
m=sippi_prior(prior);          % generate a realization from the prior model
sippi_plot_model(prior,m,im)      % visualize the realization from the prior
print_mul('prior_example_5_2d_snesim')


return

%% 1D GAUSSIAN


% Gaussian wil normal score transformation / target distributino
im=im+1;
prior{im}.name='Gaussian w. target distribution';
prior{im}=p{im-1};

N=10000;prob_chan=0.5;
d1=randn(1,ceil(N*(1-prob_chan)))*5+1900;
d2=randn(1,ceil(N*(prob_chan)))*5+2100;
d3=randn(1,ceil(N*(prob_chan)))*5+2300;
d_target=[d1(:);d2(:)];
d_target=[d1(:);d2(:);d3(:)];
[d_nscore,o_nscore]=nscore(d_target,1,1,min(d_target),max(d_target),0);
prior{im}.o_nscore=o_nscore;

% SISIM
im=im+1;
prior{im}.type='SISIM';
prior{im}.name='Sequential Indicator Simulation';
prior{im}.x=x; % X array
prior{im}.y=y_flat; % Y array
prior{im}.m0=2000;
prior{im}.Va='200 Gau(100,90,.01)';
prior{im}.Va='200 Gau(30,90,.05)';
prior{im}.marginal_prob=[0.4 0.2 .4];


% SNESIM
im=im+1;
prior{im}.type='SNESIM';
prior{im}.x=x; % X array
prior{im}.y=y_flat; % Y array
prior{im}.S=sgems_get_par('snesim_std');
prior{im}.scaling=0.5;
prior{im}.rotation=00;



% generate a realization
p=sippi_prior_init(p);
m=sippi_prior(p);
sippi_plot_model(p,m)
