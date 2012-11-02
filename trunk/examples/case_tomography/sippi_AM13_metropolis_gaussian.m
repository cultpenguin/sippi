% sippi_AM13_metropolis_gaussian.m


%% Load the travel time data set from ARREN�S
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13';

%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
% optionally use only a subset of data
%data{id}.i_use=[10:10:length(data{id}.d_obs)];
data{id}.Ct=1; % Covaiance describing modelization error

%% SETUP PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
prior{im}.x=[-1:.15:6];
prior{im}.y=[0:.15:13];

prior{im}.cax=[.1 .18];


prior{im}.seq_gibbs.type=2;

%prior{1}.seq_gibbs.step=.01;
%prior{1}.seq_gibbs.step_min=0.0001;
%prior{1}.seq_gibbs.step_max=1;

prior=sippi_prior_init(prior);
% plot a sample of the prior model
sippi_plot_prior(prior);

% generate and plot one realization of the prior model
[m,prior]=sippi_prior(prior);
sippi_plot_model(prior,m)


%% SETUP THE FORWARD MODEL
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
%forward.type='eikonal';
forward.type='ray';forward.linear=1;

% Compute the forward response related to the realization of the prior
% model generated above
[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
% plot the forward response and compare it to the observed data
sippi_plot_data(d,data); 

[logL,L,data]=sippi_likelihood(d,data);

%% SETUP METROPOLIS
options.mcmc.nite=100000;
options.mcmc.i_plot=2000;
%options.mcmc.i_sample=250;
randn('seed',2);rand('seed',2);
options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

