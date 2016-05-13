clear all;close all

forward.forward_function='sippi_forward_fmm';
forward.fmm.N=3;

% prior;
ip=1;
prior{ip}.type='fftma';

prior{ip}.x=1:1:20;
prior{ip}.y=1:1:20;
prior{ip}.m0=0;
prior{ip}.Cm='1 Gau(1,90,.2)';
prior{ip}.d_target=[0 0 0 1 1 1];
[m_ref,prior]=sippi_prior(prior);

% REFERENCE DATA
use_ti=1;
if use_ti==1;
    TI=channels;
    %TI=TI(3:3:end,3:3:end)
    m_ref{1}=TI;
end
[d,forward]=sippi_forward_fmm(m_ref,forward);

data{1}.d_obs=d{1};
data{1}.noise_model='sippi_likelihood_fmm';


% TEST
m=sippi_prior(prior);
[d,forward]=sippi_forward(m,forward);

[logL]=sippi_likelihood(d,data)

%% MCMC
prior{1}.seq_gibbs.step=0.01;
options.mcmc.nite=1000;   % [1] : Number if iterations
options.mcmc.i_sample=100; % : Number of iterations between saving model to disk
options.mcmc.i_plot=10;  % [1]: Number of iterations between updating plots
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)
