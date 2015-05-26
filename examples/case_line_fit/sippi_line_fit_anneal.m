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
prior{im}.m_true=D.m_ref{1};

% 1st order, the gradient
im=2;
prior{im}.type='gaussian';
prior{im}.name='gradient';
prior{im}.m0=0;
prior{im}.std=4;
prior{im}.norm=80;
prior{im}.m_true=D.m_ref{2};

% 2nd order
im=3;
prior{im}.type='gaussian';
prior{im}.name='2nd';
prior{im}.m0=0;
prior{im}.std=1;
prior{im}.norm=80;
prior{im}.m_true=D.m_ref{3};


%% Setup the forward model in the 'forward' structure
forward=D.forward;

%% Set up the 'data' structure
data=D.data;

%% Perform extended Metropolis sampling 
% set some MCMC options.
options.mcmc.nite=20000;
options.mcmc.i_plot=1000;
options.txt='case_line_fit';

for ip=1:length(prior)
    prior{ip}.seq_gibbs.i_update_step_max=options.mcmc.nite;
    prior{ip}.seq_gibbs.i_update_step=100;options.mcmc.nite;
end

doAnneal=1;
if doAnneal==1
    %options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=options.mcmc.nite; %  iteration number when annealing begins
    options.mcmc.anneal.fac_begin=2; % default, noise is scaled by fac_begin at iteration i_begin
    options.mcmc.anneal.fac_end=.01; % default, noise is scaled by fac_end at iteration i_end
end

[options]=sippi_metropolis(data,prior,forward,options);

for i=1:length(prior);
    disp(sprintf('%s ref=%g, true=%g ',prior{i}.name,prior{i}.m_true,options.C{1}.m_current{i}))
end
