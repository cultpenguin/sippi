% sippi_AM13_metropolis_gaussian_covariance_inference.m
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13_covinf';

%% SETUP DATA, PRIOR and FORWARD
% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=[10:10:length(data{id}.d_obs)];
data{id}.Ct=1; % modelization error
data{id}.Ct=1+D.Ct; % modelization and static error

% SETUP PRIOR

% range - horizontal
im=1;
prior{im}.type='gaussian';
prior{im}.x=1;
prior{im}.m0=8;
prior{im}.min=0;
prior{im}.max=20;
prior{im}.std=6;
prior{im}.name='range_1';
prior{im}.seq_gibbs.step_min=0.01; % MAKE SURE SIPPI_PRIOR_INIT WORKS OK
prior{im}.seq_gibbs.step_max=1;
prior{im}.seq_gibbs.step=.8;
prior{im}.prior_master=4;
prior{im}.norm=20;

% range - vertical
im=im+1;
prior{im}=prior{im-1};
prior{im}.name='range_2';


% rotation
im=im+1;
prior{im}.type='gaussian';
prior{im}.name='ang_1';
prior{im}.m0=90;
prior{im}.std=10;
prior{im}.norm=2;


% velocity field
im=im+1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
prior{im}.x=[-1:.2:6];
prior{im}.y=[0:.2:13];



prior{im}.cax=[.1 .18];

% update i_master
i_master=im;
for i=1:(i_master-1);
    prior{i}.prior_master=i_master;
    prior{i}.seq_gibbs.type=2;prior{i}.seq_gibbs.step=1;
end


prior=sippi_prior_init(prior);

% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';


%% SETUP METROPOLIS
options.mcmc.nite=500000;
options.mcmc.i_plot=200;
options.mcmc.i_sample=250;

options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

