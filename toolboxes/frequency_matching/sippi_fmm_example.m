% sippi_fmm_example
%
%
% Example replicating the example if section 3 in 
%  "
%  Cordua et al., 2015 - Improving the Pattern Reproducibility of
%  Multiple-Point-Based Prior Models Using Frequency Matching.
%  Mathematical Geosciences 47(3).
% "
%
% See also: sippi_forward_fmm, sippi_likelihood_fmm, multinomial
%
clear all;close all;

useSnesimPrior=1;
SampleOnlyPrior=0;
%step=[.25 .5];
step=[2 2];

TI=channels;
TI=TI(3:3:end,3:3:end);dx=0.25;
%TI=TI(2:2:end,2:2:end);dx=0.1;
%TI=TI(1:1:end,1:1:end);dx=0.05;
x=0:dx:5;
y=0:dx:12;

%% SETUP PRIOR
if useSnesimPrior==1;
    ip=1;
    prior{ip}.type='mps';
    prior{ip}.method='mps_snesim';
    prior{ip}.x=x;
    prior{ip}.y=y;
    prior{ip}.ti=TI;
    options.txt='TIprior';
else
    % uniform
    ip=1;
    prior{ip}.type='fftma';
    prior{ip}.x=x;
    prior{ip}.y=y;
    prior{ip}.m0=0;
    prior{ip}.Cm='1 Sph(.01)';
    prior{ip}.d_target=[0 0 0 0 0 0 0 1 1 1];
    options.txt='UNIFprior';
end
options.txt=sprintf('%s_%gcm',options.txt,dx*100);
[m,prior]=sippi_prior(prior);
%
figure(1);
subplot(1,2,1);
imagesc([1:1:size(TI,2)].*dx,[1:1:size(TI,1)].*dx,TI);
axis image
subplot(1,2,2);
imagesc(prior{1}.x,prior{1}.y,m{1});
axis image
print_mul(options.txt)

%% forward
forward.forward_function='sippi_forward_fmm';
forward.fmm.N=5;

%% reference data
m_ref{1}=TI;
[d,forward]=sippi_forward(m_ref,forward);
data{1}.d_obs=d{1}
data{1}.noise_model='sippi_likelihood_fmm';

% set noise/prior
data{1}.nprior=1;
n_av_count_per_bin=ceil(sum(d{1})./length(d{1}));
data{1}.nprior=n_av_count_per_bin*1;


%% TEST SETUP
m=sippi_prior(prior);
[d,forward]=sippi_forward(m,forward);
[logL]=sippi_likelihood(d,data)

%% MCMC
options.mcmc.accept_all=SampleOnlyPrior;
if (options.mcmc.accept_all==1)
    options.txt=[options.txt,'_OnlyPrior'];
end
prior{1}.seq_gibbs.i_update_step_max=1;
prior{1}.seq_gibbs.type=1;
prior{1}.seq_gibbs.step=step;
options.mcmc.nite=100000;   % [1] : Number if iterations
options.mcmc.i_sample=500; % : Number of iterations between saving model to disk
options.mcmc.i_plot=100;  % [1]: Number of iterations between updating plots
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)




