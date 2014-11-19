% sippi_linefit_2nd_order: Fiting line using SIPPI
clear all;close all
rand('seed',1);randn('seed',1);

%% LOAD DATA
load('sippi_linefit_data');

%% Setting up the prior model

% the intercept
im=1;
prior{im}.type='gaussian';
prior{im}.name='intercept';
prior{im}.m0=0;
prior{im}.std=30;

% 1st order, the gradient
im=2;
prior{im}.type='gaussian';
prior{im}.name='gradient';
prior{im}.m0=0;
prior{im}.std=4;
prior{im}.norm=80;

% 2nd order
im=3;
prior{im}.type='gaussian';
prior{im}.name='2nd';
prior{im}.m0=0;
prior{im}.std=1;
prior{im}.norm=80;


%% Perform extended Metropolis sampling 
options.plot.hardcopy_types=0; % NO HARDCOPY 
% set some MCMC options.
options.mcmc.nite=40000;
options.mcmc.i_sample=50;
options.mcmc.i_plot=100;
options.mcmc.m_ref=m_ref;
options.txt='case_line_fit_2nd_order';

[options]=sippi_metropolis(data,prior,forward,options);
sippi_plot_prior(options.txt);
sippi_plot_posterior(options.txt);
