% nmo_metropolis
%
% see also nmo_setup_example
%
clear all;close all;

% load data from nmo_setup_example
mat_file='nmo_reference_data_type2_SN1_nx01.mat';
if exist(mat_file);
  load(mat_file)
else
  SN=1;
  ptype=2;
  nx=1;
  nmo_setup_example;
end

try;clear options;end
for ip=1:length(prior);
  prior{ip}.seq_gibbs.i_update_step_max=10000;
end

options.mcmc.nite=10000;   % [1] : Number if iterations
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

%%
m1=sippi_get_sample(options.txt,1);
m2=sippi_get_sample(options.txt,2);
m3=sippi_get_sample(options.txt,3);
%%
for i=1:size(m1,3);
  m{1}=m1(:,i);
  m{2}=m2(:,i);
  m{3}=m3(:,i);
  d=sippi_forward(m,forward,prior);
  sippi_plot_data_reflection_nmo(d,data,1,prior);
  drawnow;
end
sippi_plot_data_reflection_nmo(d_ref,data,1,prior);
drawnow;
  
