% sippi_AM13_metropolis_bimodal 2D inversion using the extended Metropolis sampler (Bimodal prior) 
%
% Example of inverting 2D Arren�s tomographic data (AM13)
% using the extended Metropolis sampler and a 
% Bimodal a priori model
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%

clear all;close all
D=load('AM13_data.mat');
options.txt='AM13_bimodal';

%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=[10:10:length(data{id}.d_obs)];
data{id}.Ct=1; % modelization error
data{id}.Ct=1+D.Ct; % modelization and static error

%% SETUP PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
prior{im}.x=[-1:.2:6];
prior{im}.y=[0:.2:13];

prior{im}.cax=[.1 .18];

% bimodal distribution
N=10000;
prob_chan=0.5;
dd=.014*2;
d1=randn(1,ceil(N*(1-prob_chan)))*.01+0.145-dd;  %0.1125;
d2=randn(1,ceil(N*(prob_chan)))*.01+0.145+dd; %0.155;
d=[d1(:);d2(:)];
[d_nscore,o_nscore]=nscore(d,1,1,min(d),max(d),0);
prior{im}.o_nscore=o_nscore;

prior=sippi_prior_init(prior);

%% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
%forward.type='fat';
%forward.linear=1;
forward.freq=0.1;
forward.forward_function='sippi_forward_traveltime';


%% SETUP METROPOLIS
options.mcmc.nite=5000000;
options.mcmc.i_plot=5000;
options.mcmc.i_sample=2000;

options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

