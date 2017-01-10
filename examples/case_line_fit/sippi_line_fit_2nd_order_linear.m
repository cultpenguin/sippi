% sippi_linefit_2nd_order_linear: Fiting line using SIPPI / Metropolis and
% least squares
clear all;close all
rng('default')
rng(1)

%% LOAD DATA
load('sippi_linefit_data');

%% Setting up the prior model

% the intercept
im=1;
prior{im}.type='Cholesky';
prior{im}.name='intercept';
prior{im}.x=[1 2 3];
prior{im}.m0=[0];
std=[20 4 1];
prior{im}.Cmat=diag(std.^2);
m=sippi_prior(prior);

%% forward
forward.forward_function='sippi_forward_linear';
forward.G(:,1)=x.^0;
forward.G(:,2)=x.^1;
forward.G(:,3)=x.^2;

%% data
data{1}.d_obs=d_obs;
data{1}.d_std=d_std;

%% Perform extended Metropolis sampling 
options.plot.hardcopy_types=0; % NO HARDCOPY 
% set some MCMC options.
options.mcmc.nite=40000;
options.mcmc.i_sample=50;
options.mcmc.i_plot=2500;
options.mcmc.m_ref=m_ref;
options.txt='case_line_fit_2nd_order_linear';
options_metro=options;
[options_metro]=sippi_metropolis(data,prior,forward,options_metro);
sippi_plot_prior(options_metro.txt);
sippi_plot_posterior(options_metro.txt);

%% LEAST SQUARES
options_lsq=options;
options_lsq.n_reals=50;

[options_lsq]=sippi_least_squares(data,prior,forward,options_lsq);

sample_1=sippi_get_sample(options_metro.txt);
sample_2=sippi_get_sample(options_lsq.txt);
figure(1);clf,
errorbar(x,d_obs,d_std);
hold on
for i=1:size(sample_1,2);
    d_1=forward.G*sample_1(:,i);
    d_2=forward.G*sample_2(:,i);
    plot(x,d_1,'r-')
    plot(x,d_2,'g-')
end
hold off
    
    
    
    


