% sippi_AM13_metropolis_gaussian.m
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

% Next few lines for selecteing box type resim
%prior{1}.seq_gibbs.type=1;
%prior{1}.seq_gibbs.step=10;
%prior{1}.seq_gibbs.step_min=.2;
%prior{1}.seq_gibbs.step_max=10;

prior{1}.seq_gibbs.step=.01;
prior{1}.seq_gibbs.step_min=0.0001;
prior{1}.seq_gibbs.step_max=1;

prior=sippi_prior_init(prior);

% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';


%% SETUP METROPOLIS
options.mcmc.nite=500;
options.mcmc.i_plot=10;
options.mcmc.i_sample=250;
randn('seed',2);rand('seed',2);
options=sippi_metropolis(data,prior,forward,options);
return
%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

