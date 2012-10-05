% sippi_AM13_metropolis_training_image
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
im=1;
prior{im}.type='SNESIM';
prior{im}.ti='snesim_std.ti';
prior{im}.index_values=[0 1];
prior{im}.m_values=[.12 .16];
prior{im}.scaling=.7;
prior{im}.rotation=30;
prior{im}.x=[-1:.2:6];
prior{im}.y=[0:.2:13];
prior{im}.cax=[.1 .18];

prior{1}.seq_gibbs.type=1;
prior{1}.seq_gibbs.step=5;
prior{1}.seq_gibbs.step_min=0.4;
prior{1}.seq_gibbs.step_max=5;

prior=sippi_prior_init(prior);

% SETUP THE FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';


%% SETUP METROPOLIS
options.mcmc.nite=500000;
options.mcmc.i_plot=200;
options.mcmc.i_sample=250;

options=sippi_metropolis(data,prior,forward,options);

%% PLOT SAMPLE FROM PRIOR
sippi_plot_prior(options.txt);
%sippi_plot_prior(prior);

%% PLOT SAMPLE FROM POSTERIOR
sippi_plot_posterior(options.txt);

%% PLOT TI
close all

O=sgems_read(prior{1}.ti);
figure;
d_ti=O.D';
d_ti(find(d_ti==0))=prior{1}.m_values(1);
d_ti(find(d_ti==1))=prior{1}.m_values(2);
imagesc(d_ti);
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
