% sippi_linefit: Fiting line using SIPPI
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

%% TEST SETUP
m=sippi_prior(prior);
[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
[logL]=sippi_likelihood(d,data);

%% Perform extended Metropolis sampling 
% set some MCMC options.
options.mcmc.nite=40000;
%options.mcmc.T=30;
options.mcmc.adaptive_rejection=1;
options.mcmc.m_ref=m_ref;
options.txt='case_line_fit_1ord';

[options]=sippi_rejection(data,prior,forward,options);
sippi_plot_prior_sample(options.txt);
sippi_plot_posterior(options.txt);
