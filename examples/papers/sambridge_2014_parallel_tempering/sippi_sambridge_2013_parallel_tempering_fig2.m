% sippi_sambridge_2014_parallel_tempering_fig2 
% SIPPI implementation of figure 2 in 
% Sambridge, 2014. A Parallel Tempering algorithm for probabilistic
% sampling and multimodal optimizationexample from Sambridge (2013). 
% doi: 10.1093/gji/ggt342
clear all;close all

%% Define the prior
prior{1}.name='Position';
prior{1}.type='uniform';
prior{1}.min=0;
prior{1}.max=100;

%% FORWARD
forward.forward_function='sippi_forward_linear';


%% DATA
% the noise model is selected according to Figure 2
data{1}.noise_model='sippi_likelihood_bimodal_sambridge_2014';
data{1}.i_use=1;
%%
m=sippi_prior(prior);
d=sippi_forward(m,forward,prior);
logL=sippi_likelihood(d,data);

%% TRUE PD
dh=1;
hx=-10:dh:110;

dh_true=.1;
hx_true=-10:dh_true:110;
for ihx=1:length(hx_true);
    d{1}=hx_true(ihx);
    logL=sippi_likelihood(d,data);
    pd_true(ihx)=exp(logL);
end
pd_area_true=sum(pd_true)*dh_true;
pd_true=pd_true./pd_area_true;

figure(1);
bh=bar(hx_true,pd_true);
set(bh,'edgecolor','flat');
title('Reference distribution')
xlabel('Position X');
ylabel('PDF')


%%
%    pd=2.^(-hx)+2.^(-(100-hx));
%pd_sum=sum(pd);


%% EXAMPLE REJECTIONS
clear options_rejection;
options_rejection.mcmc.i_plot=1000;
options_rejection.mcmc.nite=30000;     % maximum number of iterations
options_rejection=sippi_rejection(data,prior,forward,options_rejection);
[r,em,ev,reals]=sippi_get_sample(options_rejection.txt)
%%
pd=hist(reals,hx);
pd_area=sum(pd)*dh;
pd=pd/pd_area;

figure(2);clf;
bh=bar(hx,pd);
hold on
p=plot(hx_true,pd_true,'k:','LineWidth',1)
hold off
set(bh,'edgecolor','flat','facecolor','b','FaceAlpha',.9);
title('Results from rejections sampling')
xlabel('Position X');
ylabel('PDF')
print_mul(sprintf('%s_rejection',mfilename));

%% EXAMPLE METROPOLIS
clear options_metro;
prior{1}.seq_gibbs.i_update_step_max=1;
prior{1}.seq_gibbs.step=0.2;
options_metro.mcmc.i_sample=1;
options_metro.mcmc.i_plot=1000;
options_metro.mcmc.nite=10000;     % maximum number of iterations
options_metro=sippi_metropolis(data,prior,forward,options_metro);

% get data
[reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(options_metro.txt);
[data,prior,options_metro_out,mcmc_metro]=sippi_get_posterior_data(options_metro);

%% plot
pd=hist(reals_all,hx);
pd_area=sum(pd)*dh;
pd=pd/pd_area;

figure(3);clf;

subplot(1,2,1);
bh=bar(hx,pd);
hold on
p=plot(hx_true,pd_true,'k:','LineWidth',1)
hold off
set(bh,'edgecolor','flat','facecolor','b','FaceAlpha',.9);
title(sprintf('Metropolis sampling, step=%g',prior{1}.seq_gibbs.step))
xlabel('Position X');
ylabel('PDF')

subplot(1,2,2);
plot(reals_all)
ylabel('Position X');
xlabel('Chain Step')
print_mul(sprintf('%s_metropolis',mfilename));

%% EXAMPLE METROPOLIS with PARALLEL TEMPERING
clear options_metro2;
prior{1}.seq_gibbs.i_update_step_max=1;
prior{1}.seq_gibbs.step=0.2;
options_metro2.mcmc.i_sample=1;
options_metro2.mcmc.i_plot=1000;
options_metro2.mcmc.nite=10000;     % maximum number of iterations

options_metro2.mcmc.n_chains=5; % set number of chains (def=1)
options_metro2.mcmc.T=1+[0:(options_metro2.mcmc.n_chains-1)]*.5;      % set temperature of chains [1:n_chains]
options_metro2.mcmc.chain_frequency_jump=0.1; % probability allowing a jump
                                            %  between two chains

options_metro2=sippi_metropolis(data,prior,forward,options_metro2);

% get data
[reals2,etype_mean2,etype_var2,reals_all2,reals_ite2]=sippi_get_sample(options_metro2.txt);
[data,prior,options_metro_out2,mcmc_metro2]=sippi_get_posterior_data(options_metro2);


% plot
pd=hist(reals_all2,hx);
pd_area=sum(pd)*dh;
pd=pd/pd_area;

figure(4);clf;

subplot(1,2,1);
bh=bar(hx,pd);
hold on
p=plot(hx_true,pd_true,'k:','LineWidth',1)
hold off
set(bh,'edgecolor','flat','facecolor','b','FaceAlpha',.9);
title(sprintf('Metropolis sampling w. parallel tempering, step=%g',prior{1}.seq_gibbs.step))
xlabel('Position X');
ylabel('PDF')

subplot(1,2,2);
plot(reals_all2)
ylabel('Position X');
xlabel('Chain Step')
print_mul(sprintf('%s_metropolis_parallel_temp',mfilename));
