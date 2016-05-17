clear all;close all

forward.forward_function='sippi_forward_fmm';

% prior;
ip=1;
prior{ip}.type='fftma';

prior{ip}.x=1:1:125;
prior{ip}.y=1:1:125;
%prior{ip}.x=1:1:83;
%prior{ip}.y=1:1:83;
%prior{ip}.x=1:1:33;
%prior{ip}.y=1:1:33;
prior{ip}.m0=0;
prior{ip}.Cm='1 Sph(10,90,0.001)'; % horizontal
%prior{ip}.Cm='1 Sph(.5,90,0.2)';
prior{ip}.d_target=[0 0 0 0 0 0 0 1 1 1];
[m_ref,prior]=sippi_prior(prior);

define_template=1;
if define_template==1
    T=[0 0 0;
    0 1 0;
    0 -1 0;
    0 2 0;
    0 -2 0;
    0 3 0;
    0 -3 0];
    forward.fmm.T=T;
elseif define_template==2
  T=[0 0 0;
    1 0 0;
    -1 0 0;
    2 0 0;
    -2 0 0;
    3 0 0;
    -3 0 0];
    forward.fmm.T=T;
elseif define_template==3
  T=[0 0 0;
    0 1 0;
    0 -1 0;
    0 2 0;
    0 -2 0;
    0 3 0;
    -3 0 0;
    1 0 0;
    -1 0 0;
    2 0 0;
    -2 0 0;
    3 0 0;
    -3 0 0];
    forward.fmm.T=T;
elseif define_template==4
    T=[0 0 0];
    forward.fmm.T=T;
else
    % make default template
    %forward.fmm.N=9;
    forward.fmm.N=5;
end


%% REFERENCE DATA
use_ti=1;
if use_ti==1;
    TI=channels;
    TI=TI(2:2:end,2:2:end);
    m_ref{1}=TI;
end
tic
profile on
[d,forward]=sippi_forward_fmm(m_ref,forward);
profile report
profile off
toc
data{1}.d_obs=d{1};
data{1}.noise_model='sippi_likelihood_fmm';

%% TEST
m=sippi_prior(prior);
[d,forward]=sippi_forward(m,forward);

n_av_count_per_bin=ceil(sum(d{1})./length(d{1}));
data{1}.nprior=1;
%data{1}.nprior=n_av_count_per_bin*1;
[logL]=sippi_likelihood(d,data);

%% MCMC
%prior{1}.seq_gibbs.step=0.01;
options.mcmc.nite=1000000;   % [1] : Number if iterations
options.mcmc.i_sample=5000; % : Number of iterations between saving model to disk
options.mcmc.i_plot=250;  % [1]: Number of iterations between updating plots
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)
