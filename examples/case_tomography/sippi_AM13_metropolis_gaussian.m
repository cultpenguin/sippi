% sippi_AM13_metropolis_modeling_error 2D inversion using the extended Metropolis sampler (Gaussian prior) 
%
% Example of inverting 2D Arrenaes tomographic data (AM13)
% using the extended Metropolis sampler and a 
% Gaussian a priori model
% and an approximation to the modeling error
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001


%% Load the travel time data set from ARRENAES
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13_gaussian';


%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
data{id}.Ct=1; % Covariance describing modelization error
data{id}.Ct=D.Ct; % Correlated noise model accroding to Cordua et al (2008; 2009)

%% SETUP PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
dx=0.15;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[.1 .18];


%% SETUP THE FORWARD MODEL(S)
% SETUP THE FORWARD MODEL USED IN INVERSION
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
forward.type='fat';forward.linear=1;forward.freq=0.1;


%% TEST THE SETUP 
% generate a realization from the prior
m=sippi_prior(prior);
% Compute the forward response related to the realization of the prior model generated above
[d]=sippi_forward(m,forward,prior,data);
% Compute the likelihood 
[logL,L,data]=sippi_likelihood(d,data);
% plot the forward response and compare it to the observed data
sippi_plot_data(d,data);

[logL,L,data]=sippi_likelihood(d,data);
%% SETUP METROPOLIS
options.mcmc.nite=1000000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=500;

options.mcmc.nite=100000;
options.mcmc.i_plot=10000;
options.mcmc.i_sample=100;
randn('seed',2);rand('seed',2);
options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior_sample(options.txt);

%% PLOT SAMPLE AND STATS FROM POSTERIOR
sippi_plot_posterior(options.txt);

%% PLOT PRIOR AND POSTERIO MOVIE 
sippi_plot_movie(options.txt)
