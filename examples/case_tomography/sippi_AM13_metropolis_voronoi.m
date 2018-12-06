% sippi_AM13_metropolis_vornoi 2D inversion using the extended Metropolis sampler (Gaussian prior) and a prior based on Voronoi Cells 
%
% Example of inverting 2D Arrenaes tomographic data (AM13)
% using the extended Metropolis sampler and a 
% prior based on N voronoi cells (mapped into a 2D grid using linear
% interpolation)
%
% As of 12/2014, the use of a 'voronoi' type prior model i undocumented in
% SIPPI.
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001


%% Load the travel time data set from ARRENAES
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13_voronoi';

rng('default')
rng(3);
%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
data{id}.Ct=1; % Covariance describing modelization error
data{id}.Ct=D.Ct; % Correlated noise model accroding to Cordua et al (2008; 2009)

%% SETUP PRIOR
ip=0;
cells_N_min=3;
cells_N_max=100;

ip=ip+1;
prior{ip}.type='voronoi';
dx=0.15;0.15;
prior{ip}.x=[-1:dx:6];
prior{ip}.y=[0:dx:13];
prior{ip}.cells_N=cells_N_max;
prior{ip}.cax=[.1 .18];
prior{ip}.m0=0.145;

ip=ip+1;
prior{ip}.type='uniform';
prior{ip}.name='cells_x';
prior{ip}.x=[1:cells_N_max];
prior{ip}.min=min(prior{1}.x);
prior{ip}.max=max(prior{1}.x);
prior{ip}.cax=[prior{ip}.min prior{ip}.max];
prior{ip}.prior_master=1;

ip=ip+1;
prior{ip}.type='uniform';
prior{ip}.name='cells_y';
prior{ip}.x=[1:cells_N_max];
prior{ip}.min=min(prior{1}.y);
prior{ip}.max=max(prior{1}.y);
prior{ip}.cax=[prior{ip}.min prior{ip}.max];
prior{ip}.prior_master=1;

ip=ip+1;
prior{ip}.type='fftma';
prior{ip}.name='cells_value';
prior{ip}.x=[1:cells_N_max];
prior{ip}.m0=prior{1}.m0;
prior{ip}.Va='.0003 Sph(.01)';
prior{ip}.cax=[0.1 0.18];
prior{ip}.prior_master=1;


ip=ip+1;
prior{ip}.type='uniform';
prior{ip}.name='cells_N';
prior{ip}.min=cells_N_min;
prior{ip}.max=cells_N_max;
prior{ip}.prior_master=1;



%% SETUP THE FORWARD MODEL(S)
% SETUP THE FORWARD MODEL USED IN INVERSION
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
%forward.type='fat';forward.linear=1;forward.freq=0.1;
forward.type='eikonal';
 forward.type='ray_2d';

%% TEST THE SETUP 
% generate a realization from the prior
[m,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m);figure(100);

% Compute the forward response related to the realization of the prior model generated above
[d]=sippi_forward(m,forward,prior,data);
% Compute the likelihood 
[logL,L,data]=sippi_likelihood(d,data);
% plot the forward response and compare it to the observed data
sippi_plot_data(d,data);

[logL,L,data]=sippi_likelihood(d,data);

%% SETUP METROPOLIS
options.mcmc.nite=1000000;
options.mcmc.i_plot=5000;
options.mcmc.i_sample=500;

options.mcmc.nite=150000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=100;

%options.mcmc.pert_strategy.perturb_all=1;

np=length(prior);
options.mcmc.pert_strategy.i_pert=[2:np];
options.mcmc.pert_strategy.i_pert_freq=[ones(1,np-1)./(np-1)];
%
doAnneal=0;
if doAnneal==1
    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=1000; %  iteration number when annealing stops
    options.mcmc.anneal.T_begin=20; % Starting temperature T_begin at iteration i_begin
    options.mcmc.anneal.T_end=1; % End temperature at iteration i_end
end

for ip=1:length(prior);
    prior{ip}.seq_gibbs.i_update_step_max=10000;
    prior{ip}.seq_gibbs.i_update_step=200;
end
options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior_sample(options.txt);

%% PLOT SAMPLE AND STATS FROM POSTERIOR
sippi_plot_posterior(options.txt);

return
%% SAMPLE PRIOR
prior{2}.seq_gibbs.step=0.05;
prior{3}.seq_gibbs.step=0.05;
prior{4}.seq_gibbs.step=0.001;
prior{5}.seq_gibbs.step=0.1;

for i=1:length(prior);prior{i}.perturb=0;end
prior{5}.perturb=1;
[m,prior]=sippi_prior(prior);
for i=1:1000;
    [m,prior]=sippi_prior(prior,m);
    imagesc(m{1});
    axis image;
    caxis([.11 .16]);
    drawnow;
end
