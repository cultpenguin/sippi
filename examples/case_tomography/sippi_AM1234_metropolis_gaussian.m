% sippi_AM1234_metropolis_gaussian 3D inversion using the extended Metropolis sampler (Gaussian prior) 
%
% Example of inverting 2D Arren�s tomographic data (AM13)
% using the extended Metropolis sampler and a 
% Gaussian a priori model
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%

clear all;close all
D=load('AM1234_data.mat');
options.txt='AM1234';

%% SETUP DATA, PRIOR and FORWARD
% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=[10:10:length(data{id}.d_obs)];
Ct=calc_Cd([D.S D.R],0,2^2,4^2);

data{id}.Ct=1+Ct; % modelization error

% SETUP PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
prior{im}.x=[-1:.5:6];
prior{im}.y=[-1:.5:6];
prior{im}.z=[0:.5:13];

prior{im}.cax=[.1 .18];
prior=sippi_prior_init(prior);

% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
forward.type='fat';
forward.forward_function='sippi_forward_traveltime';


m=sippi_prior(prior);
[d,forward,prior,data]=sippi_forward_traveltime(m,forward,prior,data);

%% SETUP METROPOLIS
options.mcmc.nite=2000000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=250;

options.mcmc.nite=2000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=50;


%options=sippi_metropolis_chains(data,prior,forward,options);%
%options=sippi_metropolis(data,prior,forward,options);
return
%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

