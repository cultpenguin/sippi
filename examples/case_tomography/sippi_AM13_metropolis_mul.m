% sippi_AM13_metropolis_vornoi 2D inversion using the extended Metropolis sampler (Gaussian prior) and a prior based on Voronoi Cells 
%
% Example of inverting 2D Arrenaes tomographic data (AM13)
% using the extended Metropolis sampler and a 
% prior based on N voronoi cells (mapped into a 2D grid using linear
% interpolation)
%
% As of 12/2014, the use of a 'voronoi' type prior model i undocumented in
% SIPPI.
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001


%% Load the travel time data set from ARRENAES
clear all;close all
D=load('AM13_data.mat');
options.txt='AM13';

rng('default')
rng(3);
%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
data{id}.Ct=1; % Covariance describing modelization error
data{id}.Ct=D.Ct; % Correlated noise model accroding to Cordua et al (2008; 2009)

%% SETUP MULTIPLE PRIOR MODELS PRIOR
dx=.15;
x=[-1:dx:6];
y=[0:dx:13];
cax=[0.12 0.18];
m0=0.145;

imp=0;


%% GAUSSIAN
imp=imp+1;
im=1;
P{imp}.prior{im}.type='FFTMA';
P{imp}.prior{im}.name='Gaussian_small';
P{imp}.prior{im}.m0=m0
P{imp}.prior{im}.Va='.0003 Sph(.5)';
P{imp}.prior{im}.x=x;
P{imp}.prior{im}.y=y;
P{imp}.prior{im}.cax=cax;

%% GAUSSIAN
imp=imp+1;
im=1;
P{imp}.prior{im}.type='FFTMA';
P{imp}.prior{im}.name='Gaussian';
P{imp}.prior{im}.m0=m0
P{imp}.prior{im}.Va='.0003 Sph(6,90,.5)';
P{imp}.prior{im}.x=x;
P{imp}.prior{im}.y=y;
P{imp}.prior{im}.cax=cax;

%% GAUSSIAN - Gau
imp=imp+1;
im=1;
P{imp}.prior{im}.type='FFTMA';
P{imp}.prior{im}.name='Gaussian - Gau';
P{imp}.prior{im}.m0=m0
P{imp}.prior{im}.Va='0.0000001 Nug(0) + .0003 Gau(4,90,.5)';
P{imp}.prior{im}.x=x;
P{imp}.prior{im}.y=y;
P{imp}.prior{im}.cax=cax;


%% MATERN
imp=imp+1;
im=0;
% velocity field
im=im+1;
P{imp}.prior{im}.type='FFTMA';
P{imp}.prior{im}.name='Gaussian (Matern)';
P{imp}.prior{im}.m0=m0;
P{imp}.prior{im}.Va='.0003 Mat(6,90,.3)';
P{imp}.prior{im}.x=x;
P{imp}.prior{im}.y=y;
P{imp}.prior{im}.cax=cax;
%P{imp}.prior{im}.fftma_options.pad_x=150;
%P{imp}.prior{im}.fftma_options.pad_y=150;
i_master=im;


% NU / Matern
im=im+1;
P{imp}.prior{im}.type='uniform';
P{imp}.prior{im}.name='nu';
P{imp}.prior{im}.min=0;
P{imp}.prior{im}.max=0.9;
P{imp}.prior{im}.prior_master=i_master;


%% BIMODAL
imp=imp+1;
im=1;
P{imp}.prior{im}.type='FFTMA';
P{imp}.prior{im}.name='Velocity (m/ns)';
P{imp}.prior{im}.Va='.0003 Sph(6,90,.3)';
P{imp}.prior{im}.x=x;
P{imp}.prior{im}.y=y;
P{imp}.prior{im}.cax=cax;

% bimodal distribution
N=10000;
prob_chan=0.5;
dd=.02;
d1=randn(1,ceil(N*(1-prob_chan)))*.0025+0.145-dd;  %0.1125;
d2=randn(1,ceil(N*(prob_chan)))*.0025+0.145+dd; %0.155;
d=[d1(:);d2(:)];
[d_nscore,o_nscore]=nscore(d,1,1,min(d),max(d),0);
P{imp}.prior{im}.o_nscore=o_nscore;
%forward.linear_m=0.14; 

%% PLURIGAUSSIAN
imp=imp+1;
im=1;
P{imp}.prior{im}.type='plurigaussian';
%P{imp}.prior{im}.pg_prior{1}.Cm=' 1 Gau(5,90,.5)';
%P{imp}.prior{im}.pg_map=[0.11 .11 .13 .13 .15 .15 .17 .13 .17];
P{imp}.prior{im}.pg_prior{1}.Cm=' 1 Gau(6,90,.3)';
P{imp}.prior{im}.pg_prior{2}.Cm=' 1 Sph(5,30,.4)';
P{imp}.prior{im}.pg_map=[0.11 0.11 0.11 0.13 0.13; 0.13 0.17 0.11 0.13 0.13; 0.13 0.13 0.13 0.15 0.15];
P{imp}.prior{im}.x=x;
P{imp}.prior{im}.y=y;
P{imp}.prior{im}.cax=cax;

%% VORONOI
imp=imp+1;
clear prior
ip=0;
cells_N_min=10;
cells_N_max=100;

ip=ip+1;
P{imp}.prior{ip}.type='voronoi';
P{imp}.prior{im}.name='Voronoi';
P{imp}.prior{ip}.x=x;
P{imp}.prior{ip}.y=y;
P{imp}.prior{ip}.cells_N=cells_N_max;
P{imp}.prior{ip}.cax=cax;
P{imp}.prior{ip}.m0=m0;

ip=ip+1;
P{imp}.prior{ip}.type='uniform';
P{imp}.prior{ip}.name='cells_x';
P{imp}.prior{ip}.x=[1:cells_N_max];
P{imp}.prior{ip}.min=min(x);
P{imp}.prior{ip}.max=max(x);
P{imp}.prior{ip}.cax=[min(x) max(x)];
P{imp}.prior{ip}.prior_master=1;

ip=ip+1;
P{imp}.prior{ip}.type='uniform';
P{imp}.prior{ip}.name='cells_y';
P{imp}.prior{ip}.x=[1:cells_N_max];
P{imp}.prior{ip}.min=min(y);
P{imp}.prior{ip}.max=max(y);
P{imp}.prior{ip}.cax=[min(y) max(y)];
P{imp}.prior{ip}.prior_master=1;

ip=ip+1;
P{imp}.prior{ip}.type='fftma';
P{imp}.prior{ip}.name='cells_value';
P{imp}.prior{ip}.x=[1:cells_N_max];
P{imp}.prior{ip}.m0=m0;
P{imp}.prior{ip}.Va='.0003 Sph(.01)';
P{imp}.prior{ip}.cax=[0.1 0.18];
P{imp}.prior{ip}.prior_master=1;

ip=ip+1;
P{imp}.prior{ip}.type='uniform';
P{imp}.prior{ip}.name='cells_N';
P{imp}.prior{ip}.min=cells_N_min;
P{imp}.prior{ip}.max=cells_N_max;
P{imp}.prior{ip}.prior_master=1;

%%
%PP=P;
%clear P
%P{1}=PP{3};


%% plot all prior model
for imp=1:length(P);
    sippi_plot_prior_sample(P{imp}.prior,1,5);
    print_mul(sprintf('prior_%2d_reals',imp))
end

%% SETUP THE FORWARD MODEL(S)
% SETUP THE FORWARD MODEL USED IN INVERSION
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
%forward.type='fat';forward.linear=1;forward.freq=0.1;forward.linear_m=m0;
%forward.type='ray_2d';
forward.type='eikonal';

% %% TEST THE SETUP 
% % generate a realization from the prior
% [m,prior]=sippi_prior(prior);
% sippi_plot_prior(prior,m);figure(100);
% 
% % Compute the forward response related to the realization of the prior model generated above
% [d]=sippi_forward(m,forward,prior,data);
% % Compute the likelihood 
% [logL,L,data]=sippi_likelihood(d,data);
% % plot the forward response and compare it to the observed data
% sippi_plot_data(d,data);
% 
% [logL,L,data]=sippi_likelihood(d,data);

%% SETUP METROPOLIS
n_sample=200;
options.mcmc.nite=1000000;
options.mcmc.i_plot=5000;

options.mcmc.nite=100000;
options.mcmc.i_plot=10000;

options.mcmc.i_sample=ceil(options.mcmc.nite/n_sample);

%options.mcmc.nite=100000;
%options.mcmc.i_plot=100;
%options.mcmc.i_sample=1000;

%options.mcmc.pert_strategy.perturb_all=1;

%np=length(prior);
%options.mcmc.pert_strategy.i_pert=[2:np];
%options.mcmc.pert_strategy.i_pert_freq=[ones(1,np-1)./(np-1)];
%

doTempering=1;
if doTempering==1
    %% TEMPERING
    options.mcmc.n_chains=3; % set number of chains (def=1)
    options.mcmc.T=[1 1.1 1.2];      % set temperature of chains [1:n_chains]
    options.mcmc.chain_frequency_jump=0.1; % probability allowing a jump
                                            %  between two chains
end
doAnneal=1;
if doAnneal==1
    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=ceil(options.mcmc.nite/20); %  iteration number when annealing stops
    options.mcmc.anneal.T_begin=20; % Starting temperature T_begin at iteration i_begin
    options.mcmc.anneal.T_end=1; % End temperature at iteration i_end
end

for imp=1:length(P);
    for ip=1:length(P{imp}.prior);
        P{imp}.prior{ip}.seq_gibbs.i_update_step_max=ceil(options.mcmc.nite/10);GOLF
        P{imp}.prior{ip}.seq_gibbs.i_update_step=200;
    end
end

for imp=1:length(P);
    O{imp}=options;
    O{imp}.txt=['AM13_P',num2str(imp)];
end

%%
parfor imp=1:length(P);
%for imp=1:length(P);
    O{imp}=sippi_metropolis(data,P{imp}.prior,forward,O{imp});
    sippi_plot_prior_sample(O{imp}.txt);
    sippi_plot_posterior(O{imp}.txt);    
end
save('AM13_mul')

