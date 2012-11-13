% sippi_AM13_metropolis_gaussian_covariance_inference 2D inversion with uncertain covariance model 
%
% Example of inverting 2D Arrenæs tomographic data (AM13)
% with uncertainty covariance properties. 
%
% 
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13_covinf';

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
data{id}.i_use=[10:10:length(data{id}.d_obs)];
data{id}.Ct=1; % modelization error
%data{id}.Ct=1+D.Ct; % modelization and static error

%% SETUP PRIOR
% make that any 'master' prior is defined AFTER the prior types it needs
% any FFTMA master should be defined AFTER associated covariance properties

im=0;
% range - horizontal
im=im+1;
prior{im}.type='gaussian';
prior{im}.x=1;
prior{im}.m0=6;
prior{im}.min=2;
prior{im}.max=10;
prior{im}.std=4;
prior{im}.name='range_1';
prior{im}.seq_gibbs.step_min=0.01; 
prior{im}.seq_gibbs.step_min=0.0; %
prior{im}.seq_gibbs.step_max=1;
prior{im}.seq_gibbs.step=1;
prior{im}.norm=50;

% range - vertical
im=im+1;
prior{im}=prior{im-1};
prior{im}.name='range_2';


% rotation
%im=im+1;
%prior{im}.type='gaussian';
%prior{im}.name='ang_1';
%prior{im}.m0=90;
%prior{im}.std=10;
%prior{im}.norm=2;


% velocity field
im=im+1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6,90,1)';
%prior{im}.Va='.0003 Gau(2)';
dx=0.25;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];

prior{im}.cax=[.1 .18];

% update i_master, prior{1}.prior_master = 'id of FFTMA  type prior'
i_master=im;
for i=1:(i_master-1);
    prior{i}.prior_master=i_master;
end

prior=sippi_prior_init(prior);


%% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
forward.type='fat';forward.freq=0.1;forward.linear=1;
%forward.type='born';forward.freq=0.1;%forward.linear=1;
forward.forward_function='sippi_forward_traveltime';
forward.im = i_master; % 'master' prior / the velocity --> NEEDED WHEN THERE IS MORE THE ONE A PRIORI TYPES

randn('seed',1);
rand('seed',1);

%% SETUP METROPOLIS
for im=1:length(prior)
    prior{im}.seq_gibbs.n_update_history=100;
    prior{im}.seq_gibbs.i_update_step_max=3000;
end
options.mcmc.nite=150000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=250;
options.mcmc.i_sample=50;%500;

%options.mcmc.pert_strategy.i_pert=[1 2 3];
%options.mcmc.pert_strategy.i_pert_freq=[2 2 1];

%% RUN 1, every 10th data
data{1}.i_use=[10:10:702];
try;forward=rmfield(forward,'G');end
options.txt='run1';try,forward=rmfield(forward.G);end
[o1,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
sippi_plot_posterior(o1.txt);

%% RUN 2, every 2nd data
data{1}.i_use=[5:5:702];
try;forward=rmfield(forward,'G');end
options.txt='run2';
[o2,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
sippi_plot_posterior(o2.txt);

%% RUN 3, all data
data{1}.i_use=[1:1:702];
try;forward=rmfield(forward,'G');end
options.txt='run3';
[o3,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
sippi_plot_posterior(o3.txt);
