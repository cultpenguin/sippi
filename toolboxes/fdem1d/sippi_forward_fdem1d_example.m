% sippi_forward_fdem1d_example
if ~exist('dy'); dy=1;end
if ~exist('x'); x=[0];end
if ~exist('x'); x=[0:25:100];end
if ~exist('y');y=[0:dy:124];end
if ~exist('doPlot'); doPlot=0;end
if ~exist('noise_level'); noise_level=5; end
if ~exist('noise_base'); noise_base=10; end

%% Setup prior model
im=1; % GAU
prior{im}.type='FFTMA';
prior{im}.m0=1.3;
prior{im}.Va='0.4 Gau(30,0,1)';
hx=300;
hy=40;
prior{im}.Va=sprintf('0.4 Gau(%f,90,%f)',hx,hy/hx);
%prior{im}.Va=sprintf('0.4 Sph(%f,90,%f)',hx,hy/hx);
prior{im}.x=x;
prior{im}.y=y;
prior{im}.daspect=[1 1 1];
prior{im}.cax=[0.2 3];
prior{im}.cmap=(cmap_linear([1 1 0; 0 0 0; 0 1 1]));

% 1D marginal
dstd1=0.2/2;f1=0.6;pm1=1.1;
dstd2=0.3/2;f2=0.3;pm2=2;
dstd3=0.35/2;f3=0.1;pm3=2.75;
N=10000;
dd=20*2;
d1=randn(1,ceil(N*(f1))).*dstd1+pm1;  %0.1125;
d2=randn(1,ceil(N*(f2))).*dstd2+pm2; %0.155;
d3=randn(1,ceil(N*(f3))).*dstd3+pm3; %0.155;
d_target=[d1(:);d2(:);d3(:)];

prior{im}.d_target=d_target;
prior{im}.m0=0;
sippi_plot_prior(prior);

%% Foward model
forward.forward_function='sippi_forward_fdem1d';
%if ~exist('dx_f','var'),dx_f=20;end
forward.x=prior{1}.x;
forward.ds=1;
%forward.S = readSystem('hemSystem.txt');
%forward.ndata=2*length(forward.S.freq);
forward.force_one_thread=1;


%% Create reference model and data
rng(3);
m_ref=sippi_prior(prior);plot(m_ref{1})
[d_ref,forward]=sippi_forward(m_ref,forward,prior);

% add noise 
data{1}.d_ref = d_ref{1};
data{1}.d_std = sqrt(((noise_level/100)*data{1}.d_ref).^2 + noise_base^2);
data{1}.d_noise = randn(size(data{1}.d_std)).*data{1}.d_std;
data{1}.d_obs = data{1}.d_ref + data{1}.d_noise;

options.sippi_plot_data_function='sippi_plot_data_fdem1d';
%sippi_plot_data_fdem1d(d_ref,data,1,options)
sippi_plot_data(d_ref,data,1,options)

%% Setup inversion
options.mcmc.nite=30000;   % [1] : Number if iterations
options.mcmc.i_sample=100; % : Number of iterations between saving model to disk
options.mcmc.i_plot=50;  % [1]: Number of iterations between updating plots

options.mcmc.m_ref  = m_ref;

% % ANNEALING (TEMPERATURE AS A FUNCTION OF ITERATION NUMBER)
% options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
% options.mcmc.anneal.i_end=100000; %  iteration number when annealing stops
% options.mcmc.anneal.T_begin=5; % Start temperature for annealing
% options.mcmc.anneal.T_end=1; % End temperature for annealing
% % TEMPERING
% options.mcmc.n_chains=3; % set number of chains (def=1)
% options.mcmc.T=[1 1.1 1.2];      % set temperature of chains [1:n_chains]
% options.mcmc.chain_frequency_jump=0.1; % probability allowing a jump between two chains

prior{1}.seq_gibbs.i_update_step_max=max([1000 options.mcmc.nite/20]);
prior{1}.seq_gibbs.step=1;

[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);

sippi_plot_posterior(options.txt)
                                            