% nmo_metropolis
%
% see also nmo_setup_example
%
clear all;

% load data from nmo_setup_example
%load nmo_reference_data_type2_SN20_nx01.mat
%load nmo_reference_data_type2_SN20_nx11.mat
load nmo_reference_data_type2_SN1_nx01.mat
%load nmo_reference_data_type2_SN1_nx11.mat

try;clear options;end
for ip=1:length(prior);
  prior{ip}.seq_gibbs.i_update_step_max=10000;
end

options.mcmc.nite=20000;   % [1] : Number if iterations
options.mcmc.i_sample=200; % : Number of iterations between saving model to disk
options.mcmc.i_plot=500;  % [1]: Number of iterations between updating plots

doTempering=0;
if doTempering==1
  options.mcmc.n_chains=4; % set number of chains (def=1)
  options.mcmc.T=[1 1.5 2 3];      % set temperature of chains [1:n_chains]
  options.mcmc.T=linspace(1,1.3,options.mcmc.n_chains);      % set temperature of chains [1:n_chains]
end

forward.type='akirichards';
forward.type='zoeppritz';
forward.type='weak_contrast';

%%
[m,prior]=sippi_prior(prior);
[d,forward]=sippi_forward(m,forward,prior);
[logL,L,data]=sippi_likelihood(d,data);
profile off
for i=1:100,
  [m,prior]=sippi_prior(prior,m);
  [d,forward]=sippi_forward(m,forward,prior);
  [logL,L,data]=sippi_likelihood(d,data);
end
%return
%%
[options]=sippi_metropolis(data,prior,forward,options);
