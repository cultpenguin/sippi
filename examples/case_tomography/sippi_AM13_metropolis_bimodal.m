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
prior{im}.Va='.0003 Sph(6)';
dx=0.15;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];

prior{im}.cax=[.1 .18];

% bimodal distribution
N=10000;
prob_chan=0.5;
dd=.01;
d1=randn(1,ceil(N*(1-prob_chan)))*.0025+0.145-dd;  %0.1125;
d2=randn(1,ceil(N*(prob_chan)))*.0025+0.145+dd; %0.155;
d=[d1(:);d2(:)];
[d_nscore,o_nscore]=nscore(d,1,1,min(d),max(d),0);
prior{im}.o_nscore=o_nscore;

prior=sippi_prior_init(prior);
m=sippi_prior(prior);
forward.linear_m=m{1}.*0+0.14; % Needed when m0 is set to 0 as def...

%% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
forward.type='fat';
forward.linear=1;
forward.freq=0.1;
forward.forward_function='sippi_forward_traveltime';

comp_model_error=0;
if comp_model_error==1;
    
    % SETUP THE 'OPTIMAL' FORWARD MODEL
    forward_full.forward_function='sippi_forward_traveltime';
    forward_full.sources=D.S;
    forward_full.receivers=D.R;
    forward_full.type='fat';forward_full.linear=0;forward_full.freq=0.1;
    
    % COMPUTE MODELING ERROR DUE TO USE OF forward AS OPPOSED TO forward_full
    N=600;
    [Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_full,forward,prior,data,N);
    
    % ASSIGN MODELING ERROR TO DATA
    for id=1:length(data);
        data{id}.dt=dt{id};
        data{id}.Ct=Ct{id};
    end
end



    
% TEMPERING
doTempering=1;
if doTempering==1;
    options.mcmc.n_chains=4; % set number of chains (def=1)
    options.mcmc.T=[1 1.5 2 3]; % set number of chains (def=1)
end



%% SETUP METROPOLIS
options.mcmc.nite=1000000;
options.mcmc.nite=30000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=500;

options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

%% PLOT PRIOR and POSTERIOR MOVIES
sippi_plot_movie(options.txt);

