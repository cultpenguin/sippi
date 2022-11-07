% sippi_AM13_metropolis_training_ipage 2D inversion using the extended Metropolis sampler (MPS prior)
%
% Example of inverting 2D Arrenæs tomographic data (AM13)
% using the extended Metropolis sampler and prior model 
% based on a training ipage
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%

clear all;close all
D=load('AM13_data.mat');
options.txt='AM13_ti';

%% SETUP DATA, PRIOR and FORWARD
% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=[10:10:length(data{id}.d_obs)];
data{id}.Ct=1; % modelization error
data{id}.Ct=1+D.Ct; % modelization and static error

% SETUP PRIOR
ip=1;
prior{ip}.type='mps';
prior{ip}.method='mps_snesim_tree';   
prior{ip}.method='mps_genesim';   
ti=channels;
ti=ti(5:5:end,5:5:end);
prior{ip}.ti=channels;
prior{ip}.index_values=[0 1];
prior{ip}.m_values=[.14 .1575];
prior{ip}.x=[-1:.2:6];
prior{ip}.y=[0:.2:13];
prior{ip}.cax=[.12 .16];

prior{ip}.seq_gibbs.type=2;
prior{ip}.seq_gibbs.step=1;
prior{ip}.seq_gibbs.step_min=0.95;
prior{ip}.seq_gibbs.step_max=1;

prior=sippi_prior_init(prior);

% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
forward.type='ray_2d';
forward.forward_function='sippi_forward_traveltime';

%% TEST SETUP
m=sippi_prior(prior);
d=sippi_forward(m,forward,prior);
sippi_plot_data(d,data)

%% SETUP METROPOLIS
options.mcmc.nite=50000;
options.mcmc.i_plot=50;
options.mcmc.i_sample=250;

 % ANNEALING (TEMPERATURE AS A FUNCTION OF ITERATION NUMBER)
 options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
 options.mcmc.anneal.i_end=1000; %  iteration number when annealing stops
 options.mcmc.anneal.T_begin=20; % Start temperature for annealing
 options.mcmc.anneal.T_end=1; % End temperature for annealing
 

options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

%% PLOT TI
close all

figure;
imagesc(ti);
xlabel('X pixel #')
ylabel('Y pixel #')
axis image
caxis(prior{1}.cax);
cb=colorbar;
set(get(cb,'Ylabel'),'String','Velocity (m/ns)');
set(get(cb,'Ylabel'),'FontSize',14);
ppp(10,10,12,4,1)
%colorbar_shift
print_mul('figure_ti_ref')
