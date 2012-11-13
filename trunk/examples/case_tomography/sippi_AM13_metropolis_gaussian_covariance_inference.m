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


% velocity field
im=im+1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6,90,1)';
prior{im}.Va='.0003 Sph(6)';
%prior{im}.Va='.0003 Gau(2)';
dx=0.25;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[.1 .18];
i_master=im;

% range - horizontal
im=im+1;
prior{im}.type='gaussian';
prior{im}.x=1;
prior{im}.m0=6;
prior{im}.min=1.5;
prior{im}.max=10.5;
prior{im}.std=4;
prior{im}.name='range_1';
prior{im}.seq_gibbs.step_min=0.01; 
prior{im}.seq_gibbs.step_min=0.0; %
prior{im}.seq_gibbs.step_max=1;
prior{im}.seq_gibbs.step=1;
prior{im}.norm=50;
prior{im}.prior_master=i_master;

%% range - vertical
%im=im+1;
%prior{im}=prior{im-1};
%prior{im}.name='range_2';


% rotation
%im=im+1;
%prior{im}.type='gaussian';
%prior{im}.name='ang_1';
%prior{im}.m0=90;
%prior{im}.std=10;
%prior{im}.norm=2;

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
    prior{im}.seq_gibbs.n_update_history=200;
    prior{im}.seq_gibbs.i_update_step_max=4000;
    prior{im}.seq_gibbs.i_update_step=100;
end
options.mcmc.nite=20000;150000;
options.mcmc.i_plot=5000;
options.mcmc.i_sample=50;

%options.mcmc.pert_strategy.i_pert=[1 2 3];
%options.mcmc.pert_strategy.i_pert_freq=[1 1 1];

%% RUN 1, every 10th data
close all
data{1}.i_use=[20:20:702];
try;forward=rmfield(forward,'G');end
options.txt='run1';try,forward=rmfield(forward.G);end
[o1,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
sippi_plot_posterior(o1.txt);
options.mcmc.m_init=m_current;
return
 
%% RUN 2, every 2nd data
close all;
data{1}.i_use=[5:5:702];
try;forward=rmfield(forward,'G');end
options.txt='run2';
[o2,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
sippi_plot_posterior(o2.txt);
options.mcmc.m_init=m_current;

%% RUN 3, all data
options.mcmc.nite=100000;
data{1}.i_use=[1:1:702];
try;forward=rmfield(forward,'G');end
options.txt='run3';
[o3,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
sippi_plot_posterior(o3.txt);
