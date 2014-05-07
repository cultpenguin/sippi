% Fiting line using SIPPI, simulated annealing
clear all;close all
rand('seed',1);randn('seed',1);

%% LOAD DATA
D=load('sippi_linefit_data');

%% Setting up the prior model

% the intercept
im=1;
prior{im}.type='gaussian';
prior{im}.name='intercept';
prior{im}.m0=0;
prior{im}.std=30;
prior{im}.m_true=D.intercept;

% 1st order, the gradient
im=2;
prior{im}.type='gaussian';
prior{im}.name='gradient';
prior{im}.m0=0;
prior{im}.std=4;
prior{im}.norm=80;
prior{im}.m_true=D.grad;

% 2nd order
im=3;
prior{im}.type='gaussian';
prior{im}.name='2nd';
prior{im}.m0=0;
prior{im}.std=1;
prior{im}.norm=80;
prior{im}.m_true=D.poly2;


%% Setup the forward model in the 'forward' structure
forward.x=D.x;
forward.forward_function='sippi_forward_linefit';

%% Set up the 'data' structure
data{1}.d_obs=D.d_obs;
data{1}.d_std=D.d_std;

%% Perform extended Metropolis sampling 
% set some MCMC options.
options.mcmc.nite=2000;
options.mcmc.i_plot=300;
options.txt='case_line_fit';

    
%options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
options.mcmc.anneal.i_end=options.mcmc.nite; %  iteration number when annealing begins
options.mcmc.anneal.fac_begin=2; % default, noise is scaled by fac_begin at iteration i_begin
options.mcmc.anneal.fac_end=.01; % default, noise is scaled by fac_end at iteration i_end
 
 
[options]=sippi_metropolis(data,prior,forward,options);

for i=1:length(prior);
    disp(sprintf('%s ref=%g, true=%g ',prior{i}.name,prior{i}.m_true,options.mcmc.m_current{i}))
end
