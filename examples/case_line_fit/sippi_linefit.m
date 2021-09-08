% sippi_linefit: Fiting line using SIPPI
clear all;close all
rng('default')
rng(1);
%% Load or generate data data
nd=11;
% 
mat_datfile = sprintf('sippi_linefit_data_%d',nd);
if exist([mat_datfile,'.mat'])
    load(mat_datfile);
else
    sippi_linefit_make_data;
end
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
options.mcmc.nite=40000;  % Run for 40000 iterations
options.mcmc.i_sample=100; % Save every 100th visited model to disc
options.mcmc.i_plot=5000; % Plot the progress information for every 2500 iterations
options.txt='case_line_fit_2nd_order'; % descriptive name for the output folder

%% metropolis
[options_metropolis]=sippi_metropolis(data,prior,forward,options);
out_folder = options_metropolis.txt;
sippi_plot_prior_sample(out_folder);
sippi_plot_posterior(out_folder);

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



