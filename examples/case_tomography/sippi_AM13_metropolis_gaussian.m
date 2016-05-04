% sippi_AM13_gaussian 2D inversion using the extended Metropolis sampler (Gaussian prior) 
%
% Example of inverting 2D Arrenaes tomographic data (AM13)
% using the extended Metropolis sampler and a 
% Gaussian a priori model.
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
% 
% See also: sippi_metropolis
%


%% Load the travel time data set from ARRENAES
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13_gaussian';


%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std.*0+0.4;;
data{id}.i_use=1:20;
%data{id}.Ct=1; % Data covariance describing modelization error
data{id}.Ct=D.Ct; % Correlated noise model according to Cordua et al (2008; 2009)
options.txt=[options.txt,'_noCt'];

%% SETUP PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
dx=0.15;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[.04 .18];

d_target=[randn(1,100)*.003+0.11 randn(1,100)*.003+0.16];
prior{im}.d_target=d_target;
prior{im}.m0=0; %% MAKE SURE sippi_forward_traveltime tests for a non-zero velocity

sippi_plot_prior(prior);


%% SETUP THE FORWARD MODEL(S)
% SETUP THE FORWARD MODEL USED IN INVERSION
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
forward.type='fat';forward.linear=1;forward.freq=0.1;
%forward.type='eikonal';
forward.m0=0.10;


%% TEST THE SETUP 
% generate a realization from the prior
[m,prior]=sippi_prior(prior);
% Compute the forward response related to the realization of the prior model generated above
[d,forward]=sippi_forward(m,forward,prior,data);
% Compute the likelihood 
[logL,L,data]=sippi_likelihood(d,data);
% plot the forward response and compare it to the observed data
sippi_plot_data(d,data);
return
[logL,L,data]=sippi_likelihood(d,data);
%% SETUP METROPOLIS
options.mcmc.nite=100000;
options.mcmc.i_plot=500;
options.mcmc.i_sample=50;
randn('seed',2);rand('seed',2);

% ANNEALING 
% example of starting with a high temperature, that allow higher
% exploration. The temperature is lowered to T=1, after which the 
% algorithm proceeds as a usual Metropolis sampler 
doAnneal=1;
if doAnneal==1;
    i_stop_anneal=1000;
    for im=1:length(prior);
        prior{im}.seq_gibbs.i_update_step_max=2*i_stop_anneal;
    end
    options.mcmc.anneal.T_begin=10; % Start temperature
    options.mcmc.anneal.T_end=1; % End temperature, at ptions.mcmc.anneal.
    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=i_stop_anneal; %  iteration number when annealing stops
end
    
% TEMPERING
% example of using parallel tempering (Sambridge, 2013)
doTempering=1;
if doTempering==1;
    options.mcmc.n_chains=4; % set number of chains (def=1)
    options.mcmc.T=[1 1.5 2 3]; % set number of chains (def=1)
end

options=sippi_metropolis(data,prior,forward,options);
options.mcmc.time_elapsed_in_seconds

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior_sample(options.txt);

%% PLOT SAMPLE AND STATS FROM POSTERIOR
sippi_plot_posterior(options.txt);

%% PLOT PRIOR AND POSTERIO MOVIE 
sippi_plot_movie(options.txt)
