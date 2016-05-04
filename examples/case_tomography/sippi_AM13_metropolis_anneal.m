% sippi_AM13_metropolis_gaussian 2D inversion using the extended Metropolis sampler (Gaussian prior) 
%
% Example of inverting 2D Arren�s tomographic data (AM13)
% using the extended Metropolis sampler and a 
% Gaussian a priori model
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%

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
data{id}.i_use=[10:10:length(data{id}.d_obs)];
%data{id}.Ct=D.Ct; % Covaiance describing modelization error
%data{id}.Ct=D.Ct+1; % Covaiance describing modelization error

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
% sippi_plot_prior(prior);

% generate and plot one realization of the prior model
[m,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m)


%% SETUP THE FORWARD MODEL
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
%forward.type='ray';forward.linear=1;
%forward.type='fat';forward.linear=1;forward.freq=0.1;
%forward.type='born';forward.linear=1;forward.freq=0.1;

% Compute the forward response related to the realization of the prior
% model generated above

[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
% plot the forward response and compare it to the observed data
sippi_plot_data(d,data); 

[logL,L,data]=sippi_likelihood(d,data);

%
sippi_plot_data(d,data)
%% SETUP METROPOLIS/ANNEALING
options.mcmc.nite=500000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=500;
randn('seed',1);rand('seed',1);

% OPTIMIZATION
options.mcmc.i_plot=100;
options.mcmc.nite=100000;
options.mcmc.i_sample=100;
options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
options.mcmc.anneal.i_end=options.mcmc.nite; %  iteration number when annealing stops
options.mcmc.anneal.T_begin=4; % start temperature at iteration i_begin
options.mcmc.anneal.T_end=.1; % end temperature iteration i_end
for i=1:length(prior);
    prior{i}.seq_gibbs.i_update_step_max=options.mcmc.nite;
end

% ONLY ANNEALING IN AN INITIAL (BURNIN) PHASE
% options.mcmc.nite=10000;
% options.mcmc.i_plot=100;
% options.mcmc.i_sample=25;
% options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
% options.mcmc.anneal.i_end=2000; %  iteration number when annealing stops
% options.mcmc.anneal.T_begin=2; % start temperature iteration i_begin
% options.mcmc.anneal.T_end=1; % end tenperature at iteration i_end

% START SAMPLING/OPTIMIZATION
options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

%% PLOT PRIOR AND POSTERIOR MOVIE
sippi_plot_movie(options.txt,1,1000,0);

%% PLOT LAST MODEL / 'OPTIMAL' MODEL
clf;sippi_plot_prior(prior,options.mcmc.m_current);
print_mul(sprintf('%s_last_model',options.txt))

[d,forward,prior,data]=sippi_forward(options.mcmc.m_current,forward,prior,data);
sippi_plot_data(d,data); 



