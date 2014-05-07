% sippi_linefit_2nd_order: Fiting line using SIPPI
clear all;close all
rand('seed',1);randn('seed',1);

%% Setting up the prior model

% the intercept
im=1;
prior{im}.type='gaussian';
prior{im}.name='intercept';
prior{im}.m0=0;
prior{im}.std=30;
prior{im}.m_true=-30;

% 1st order, the gradient
im=2;
prior{im}.type='gaussian';
prior{im}.name='gradient';
prior{im}.m0=0;
prior{im}.std=4;
prior{im}.norm=80;
prior{im}.m_true=2;

% 2nd order
im=3;
prior{im}.type='gaussian';
prior{im}.name='2nd';
prior{im}.m0=0;
prior{im}.std=1;
prior{im}.norm=80;
prior{im}.m_true=.9;

%% Setup the forward model in the 'forward' structure
nd=40;
forward.x=linspace(1,20,nd)';
forward.forward_function='sippi_forward_linefit';

%% Set up the 'data' structure
for ip=1:length(prior);
    m_true{ip}=prior{ip}.m_true;
end
d=sippi_forward_linefit(m_true,forward);
d_obs=d{1};
% Add noise top data
data{1}.d_std=10;
data{1}.d_obs=d_obs+randn(size(d_obs)).*data{1}.d_std;

%% Perform extended Metropolis sampling 
% set some MCMC options.
options.mcmc.nite=40000;
options.mcmc.i_sample=50;
options.mcmc.i_plot=2500;
options.txt='case_line_fit_2nd_order';

[options]=sippi_metropolis(data,prior,forward,options);
sippi_plot_prior_sample(options.txt);
sippi_plot_posterior(options.txt);
