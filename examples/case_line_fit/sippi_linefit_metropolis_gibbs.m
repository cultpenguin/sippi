% sippi_linefit: Fiting line using SIPPI
clear all;close all
rng('default')
%rng(1);
%% LOAD DATA
%nd=21; sippi_linefit_make_data;
load('sippi_linefit_data_11');

figure(1);clf;
e=errorbar(x,d_obs,d_std,'k*');
set(e,'LineStyle','none')
xlabel('x')
ylabel('d')
box on
grid on
axis([min(x)-1 max(x)+1 min(d_obs)-20 max(d_obs)+20])
ppp(10,4,12,2,2)
print_mul(sprintf('sippi_linefit_data_%d',length(d_obs)));

%% setup the data
data{1}.d_obs=d_obs;
data{1}.d_std=d_std;

%% setup the forward
forward.x=x;
forward.forward_function='sippi_forward_linefit';


%% setup the prior model
% the intercept
im=1;
prior{im}.type='gaussian';
prior{im}.name='intercept';
prior{im}.m0=0;
prior{im}.std=30;
prior{im}.m_true=m_ref{1};

% 1st order, the gradient
im=2;
prior{im}.type='gaussian';
prior{im}.name='gradient';
prior{im}.m0=0;
prior{im}.std=4;
prior{im}.norm=80;
prior{im}.m_true=m_ref{2};

% 2ndst order
im=3;
prior{im}.type='uniform';
prior{im}.name='2nd';
prior{im}.min=-.3;
prior{im}.max=.3;
prior{im}.m_true=m_ref{3};


%% TEST SETUP
m=sippi_prior(prior);
[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
[logL]=sippi_likelihood(d,data);

%%
figure(1);clf;
e=errorbar(x,d_obs,d_std,'k*');
hold on
plot(forward.x,d{1},'r*')
hold off
set(e,'LineStyle','none')
xlabel('x')
ylabel('d')
box on
grid on
axis([-1 21 -60 40])
ppp(10,4,12,2,2)
print_mul(sprintf('sippi_linefit_data_%d_prior',length(d_obs)));



%% Sample the posterior
% Perform extended Metropolis sampling 
options.mcmc.nite=4000;  % Run for 40000 iterations
options.mcmc.i_sample=10; % Save every 100th visited model to disc
options.mcmc.i_plot=1000; % Plot the progress information for every 2500 iterations
options.mcmc.i_plot=10;
%options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
%options.mcmc.anneal.i_end=options.mcmc.nite; %  iteration number when annealing stops
%options.mcmc.anneal.T_begin=15; % Start temperature for annealing
%options.mcmc.anneal.T_end=.01; % End temperature for annealing
for ip=1:length(prior);prior{ip}.seq_gibbs.i_update_step_max=1000;end


%% TEMPERING
%options.mcmc.n_chains=3; % set number of chains (def=1)
%options.mcmc.T=[1 2 8];      % set temperature of chains [1:n_chains]
%options.mcmc.chain_frequency_jump=0.1; % probability allowing a jump
                                            %  between two chains

options.txt='case_line_fit_2nd_order'; % descriptive name for the output folder


%% PERTUBATION STRAT
options.mcmc.pert_strategy.i_pert = [1]; % only perturb prior 1
%options.mcmc.pert_strategy.i_pert_freq = [1]; % only perturb prior 1

options.mcmc.gibbs.usedim = 2; % use 1D or 2D gibbs sampling
options.mcmc.gibbs.i_pert = [2 3]; % select the prior ids to use for Gibbs sampling (must be 1D)
options.mcmc.gibbs.i_gibbs = 50; % perform gibbs sampling for eevery i_gibbs iterations!options.mcmc.gibbs.Nn2=11;;


options.mcmc.gibbs.Nn1=51;;
options.mcmc.gibbs.Nn2=51;;
%% metropolis
[options_metropolis]=sippi_metropolis_gibbs(data,prior,forward,options);
%[options_metropolis]=sippi_metropolis(data,prior,forward,options);
%sippi_plot_prior_sample(options_metropolis.txt);
sippi_plot_posterior_2d_marg(options_metropolis.txt);

return
%% rejection
options.mcmc.adaptive_rejection=1;
[options_rejection]=sippi_rejection(data,prior,forward,options);
sippi_plot_prior_sample(options_rejection.txt);
sippi_plot_posterior(options_rejection.txt);


%% Plot 1D histograms
for i=1:3;
    [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(options_metropolis.txt,i);
    m_post(:,i)=reals_all;
end
figure(2);clf;
e=errorbar(x,d_obs,d_std,'k*');
hold on
for i=1:size(m_post,1);
    m{1}=m_post(i,1);
    m{2}=m_post(i,2);
    m{3}=m_post(i,3);
    [d,forward,prior,data]=sippi_forward(m,forward,prior,data);
    plot(x,d{1},'k-','color',[1 1 1].*.8);
end
[d,forward,prior,data]=sippi_forward(m_ref,forward,prior,data);
plot(x,d{1},'g-','LineWidth',3);
e=errorbar(x,d_obs,d_std,'k*');
hold off

set(e,'LineStyle','none')
xlabel('x')
ylabel('d')
box on
grid on
axis([-1 21 -60 40])
ppp(10,4,12,2,2)
print_mul(sprintf('sippi_linefit_data_%d_post',length(d_obs)));



