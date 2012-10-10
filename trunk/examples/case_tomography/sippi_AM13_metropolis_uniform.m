% sippi_AM13_metropolis_uniform.m
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13_uniform';

%% SETUP DATA, PRIOR and FORWARD
% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=[1:1:length(data{id}.d_obs)];
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

% uniform distribution
N=10000;
d=(rand(1,N)-.5)*.09+.145;
[d_nscore,o_nscore]=nscore(d,1,1,min(d),max(d),0);
prior{im}.o_nscore=o_nscore;

prior=sippi_prior_init(prior);

% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
forward.forward_function='sippi_forward_traveltime';


%% SETUP METROPOLIS
options.mcmc.nite=500000;
options.mcmc.i_plot=200;
options.mcmc.i_sample=250;

options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

