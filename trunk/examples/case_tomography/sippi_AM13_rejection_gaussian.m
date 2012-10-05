% sippi_AM13_rejection_gaussian.m
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13';

%% SETUP DATA, PRIOR and FORWARD
% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=[10:10:length(data{id}.d_obs)];
data{id}.Ct=1; % modelization error
data{id}.Ct=1+D.Ct; % modelization and static error

% SETUP PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
prior{im}.x=[-1:.2:6];
prior{im}.y=[0:.2:13];

prior{im}.cax=[.1 .18];
prior=sippi_prior_init(prior);

% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';


%% SETUP METROPOLIS
options.mcmc.nite=500000;
options.mcmc.i_plot=200;
options.mcmc.i_sample=250;

options=sippi_rejection(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

