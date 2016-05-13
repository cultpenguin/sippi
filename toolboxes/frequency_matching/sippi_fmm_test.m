clear all;close all

forward.forward_function='sippi_forward_fmm';
forward.fmm.N=3;

% REFERENCE DATA

TI=channels;
TI=TI(4:4:end,4:4:end)
m_ref{1}=TI;
[d,forward]=sippi_forward_fmm(m_ref,forward);

data{1}.d_obs=d{1};
data{1}.noise_model='sippi_likelihood_fmm';
% prior;
ip=1;
prior{ip}.type='fftma';


prior{ip}.x=1:1:20;
prior{ip}.y=1:1:20;
prior{ip}.m0=0;
prior{ip}.Cm='1 Sph(.1,90,.2)';
prior{ip}.d_target=[0 0 0 1 1 1];
[m,prior]=sippi_prior(prior);

% TEST
m=sippi_prior(prior);
[d,forward]=sippi_forward(m,forward);

[logL]=sippi_likelihood(d,data);

%% MCMC
options.mcmc.nite=1000;   % [1] : Number if iterations
options.mcmc.i_sample=100; % : Number of iterations between saving model to disk
options.mcmc.i_plot=10;  % [1]: Number of iterations between updating plots
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)
