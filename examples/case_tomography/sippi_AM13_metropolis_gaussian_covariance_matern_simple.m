% sippi_AM13_metropolis_gaussian_covariance_matern_simple 2D inversion with uncertain covariance model 
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
data{id}.Ct=1; % modelization error
%data{id}.Ct=1+D.Ct; % modelization and static error

%% SETUP PRIOR
im=0;

% velocity field
im=im+1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Mat(6,90,.3)';
dx=0.25;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[.1 .18];
prior{im}.fftma_options.pad_x=150;
prior{im}.fftma_options.pad_y=150;
i_master=im;


% NU / Matern
im=im+1;
prior{im}.type='uniform';
prior{im}.name='nu';
prior{im}.min=0;
prior{im}.max=0.9;
prior{im}.prior_master=i_master;

prior=sippi_prior_init(prior);

%sippi_plot_prior_sample(prior,1);
%% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
%forward.type='eikonal';
forward.type='fat';forward.freq=0.1;forward.linear=1;
%forward.type='born';forward.freq=0.1;%forward.linear=1;
forward.forward_function='sippi_forward_traveltime';
forward.im = i_master; % 'master' prior / the velocity --> 
                       % NEEDED WHEN THERE IS MORE THE ONE A PRIORI TYPES
                       % (def, im=1);
rng('default')
rng(1);
                       

%% SETUP METROPOLIS
for im=1:length(prior) 
    prior{im}.seq_gibbs.n_update_history=200;% optional
    prior{im}.seq_gibbs.i_update_step_max=4000;% optional
    prior{im}.seq_gibbs.i_update_step=100;% optional
end
options.mcmc.nite=1000000;% optional
options.mcmc.i_plot=1000;% optional
options.mcmc.i_sample=250;% optional

[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
sippi_plot_posterior(options.txt);
