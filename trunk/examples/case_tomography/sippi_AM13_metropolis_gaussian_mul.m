% sippi_AM13_metropolis_gaussian 2D inversion using the extended Metropolis sampler (Gaussian prior) 
%
% Example of inverting 2D Arrenæs tomographic data (AM13)
% using the extended Metropolis sampler and a 
% Gaussian a priori model
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%

%% Load the travel time data set from ARRENÆS
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
%data{id}.Ct=D.Ct; % Covaiance describing modelization error
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

prior=sippi_prior_init(prior);

% plot a sample of the prior model
sippi_plot_prior(prior);

% generate and plot one realization of the prior model
[m,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m)

%% SETUP THE FORWARD MODEL
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
%forward.type='eikonal';
%forward.type='ray';forward.linear=1;
forward.type='fat';forward.linear=1;forward.freq=0.1;
%forward.type='born';forward.linear=1;forward.freq=0.1;

% Compute the forward response related to the realization of the prior
% model generated above

[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
[logL,L,data]=sippi_likelihood(d,data);

%%
make_synth=1;
if make_synth==1;
    d_noise=gaussian_simulation_cholesky(0,data{1}.CD,1);
    data{1}.d_obs=d{1}+d_noise;
    [d,forward,prior,data]=sippi_forward(m,forward,prior,data);
end

% plot the forward response and compare it to the observed data
sippi_plot_data(d,data); 

[logL,L,data]=sippi_likelihood(d,data);

%% SETUP METROPOLIS
options.mcmc.nite=500000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=100;
randn('seed',2);rand('seed',2);
for i=1:3;
  o{i}=sippi_metropolis(data,prior,forward,options);
  % PLOT SAMPLE FROM PRIOR
  sippi_plot_prior_sample(o{i}.txt);
  % PLOT SAMPLE AND STATS FROM POSTERIOR
  sippi_plot_posterior(o{i}.txt);
end

