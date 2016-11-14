% sippi_AM13_gaussian 2D inversion using the extended Metropolis sampler with uncertain noise model (Gaussian prior) 
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
clear all;close all;
rng('default');rng(1);
D=load('AM13_data.mat');
options.txt='AM13_gaussian_noise';
%D2=D;
%D.S(352:end,:)=D2.R(352:end,:);
%D.R(352:end,:)=D2.S(352:end,:);

%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std.*0+0.4;;
data{id}.Ct=D.Ct; % Correlated noise model according to Cordua et al (2008; 2009)
% If data{id}.Ct, then the magniatude of an uncorrelated Gaussian noise, on
% top of the modeling error is inferrred


%% SETUP PRIOR
dx=.15;
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[-1 1].*.035+prior{im}.m0;

inferNoise=1;
if inferNoise==1;
    % in order to infer the the magnitude of uncorrelated measurement noise
    % define a 1D prior with name 'd_std' (or 'd_var'), and set  prior{im}.perturb_noise=1;
    % as for example in the case below:
    im=im+1;
    prior{im}.type='uniform';
    prior{im}.name='d_std';
    prior{im}.perturb_noise=1; % 
    %prior{im}.perturb_noise_i_data=1; %def=1
    prior{im}.min=.1;
    prior{im}.max=4;
    prior{im}.seq_gibbs.step_min=0.01;   
    profile on

    options.txt=sprintf('%s_noise',options.txt);
end



%% SETUP THE FORWARD MODEL(S)
% SETUP THE FORWARD MODEL USED IN INVERSION
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
%forward.type='fat';forward.linear=1;forward.freq=0.1;
%forward.type='ray';
forward.type='ray_2d';forward.r=2;
%forward.type='eikonal';
%forward.type='fd';


%% TEST THE SETUP 
% generate a realization from the prior
[m,prior]=sippi_prior(prior);

% Compute the forward response related to the realization of the prior model generated above
[d,forward]=sippi_forward(m,forward,prior,data);
try;
    forward.G=sparse(forward.G);
end

% plot the geomoetry
sippi_plot_traveltime_kernel(forward,prior,m);
figure(1);print_mul('AM13_SourceReceiverLoc');
figure(2);print_mul('AM13_Kernel');

% Compute the likelihood 
[logL,L,data]=sippi_likelihood(d,data);
% plot the forward response and compare it to the observed data
sippi_plot_data(d,data);
print_mul('AM13_data');


%% SETUP METROPOLIS
options.mcmc.nite=500000;
options.mcmc.i_plot=10000;
n_reals_out=200;
options.mcmc.i_sample=options.mcmc.nite/n_reals_out;
rng(1);

% ANNEALING 
% example of starting with a high temperature, that allow higher
% exploration. The temperature is lowered to T=1, after which the 
% algorithm proceeds as a usual Metropolis sampler 
doAnneal=1;
if doAnneal==1;
    i_stop_anneal=4000;
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
doTempering=0;
if doTempering==1;
    options.mcmc.n_chains=4; % set number of chains (def=1)
    options.mcmc.T=[1 1.5 2 3]; % set number of chains (def=1)
end

options.mcmc.pert_strategy.i_pert = [1,2]; % Perturb both prior 1 and 2
options.mcmc.pert_strategy.i_pert_freq = [9 1]; % Perturb prior 1 90% of the time



options=sippi_metropolis(data,prior,forward,options);

options.mcmc.time_elapsed_in_seconds

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior_sample(options.txt);

%% PLOT SAMPLE AND STATS FROM POSTERIOR
sippi_plot_posterior(options.txt);

%% PLOT PRIOR AND POSTERIO MOVIE 
%sippi_plot_movie(options.txt)
