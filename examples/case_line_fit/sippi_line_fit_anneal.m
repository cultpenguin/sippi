% Fiting line using SIPPI, simulated annealing
clear all;close all
rand('seed',1);randn('seed',1);

%% Setting up the prior model

% the gradient
im=1;
prior{im}.type='gaussian';
prior{im}.name='gradient';
prior{im}.m0=0;
prior{im}.std=4;
prior{im}.norm=80;
prior{im}.m_true=2;

% the intercept
im=im+1;
prior{im}.type='gaussian';
prior{im}.name='intercept';
prior{im}.m0=0;
prior{im}.std=30;
prior{im}.m_true=-30;

prior=sippi_prior_init(prior);

%% Setup the forward model in the 'forward' structure
nd=40;
forward.x=linspace(1,20,nd)';
forward.forward_function='sippi_forward_linefit';

%% Set up the 'data' structure
id=1;
m_true{1}=prior{1}.m_true;
m_true{2}=prior{2}.m_true;
d=sippi_forward_linefit(m_true,forward);
d_obs=d{1};
% Add noise top data
data{1}.d_std=10;
data{1}.d_obs=d_obs+randn(size(d_obs)).*data{1}.d_std;

%% Perform extended Metropolis sampling 
% set some MCMC options.
options.mcmc.nite=20000;
options.mcmc.i_sample=50;
options.mcmc.i_plot=1000;
options.txt='case_line_fit';

options.mcmc.nite=8000;
    
%options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
options.mcmc.anneal.i_end=options.mcmc.nite; %  iteration number when annealing begins
options.mcmc.anneal.fac_begin=2; % default, noise is scaled by fac_begin at iteration i_begin
options.mcmc.anneal.fac_end=.01; % default, noise is scaled by fac_end at iteration i_end
 
 
[options]=sippi_metropolis(data,prior,forward,options);

for i=1:length(prior);
disp(sprintf('%s ref=%g, optimal=%g ',prior{i}.name,prior{i}.m_true,options.mcmc.m_current{i}))
end
