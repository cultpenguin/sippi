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
options.txt='AM13_modeling_error';


%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
data{id}.Ct=D.Ct+1; % Covaiance describing modelization error

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
forward.type='eikonal'
forward.type='fat';forward.linear=1;forward.freq=0.1;

% SETUP THE 'OPTIMAL' FORWARD MODEL
forward_full.forward_function='sippi_forward_traveltime';
forward_full.sources=D.S;
forward_full.receivers=D.R;
forward_full.type='fat';forward_full.linear=0;forward_full.freq=0.1;

% COMPUTE MODELING ERROR DUE TO USE OF forward AS OPPOSED TO forward_full
N=600;
[Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_full,forward,prior,N);

% ASSIGN MODELING ERROR TO DATA
for id=1:length(data);
  data{id}.dt=dt{id};
  data{id}.Ct=Ct{id};
end 

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
randn('seed',2);rand('seed',2);
options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior_sample(options.txt);

%% PLOT SAMPLE AND STATS FROM POSTERIOR
sippi_plot_posterior(options.txt);


