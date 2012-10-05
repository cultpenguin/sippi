% Fiting line using SIPPI
clear all;close all
rand('seed',1);randn('seed',1);

%% Setting up the prior model

% the gradient
im=1;
prior{im}.type='gaussian';
prior{im}.name='gradient';
prior{im}.m0=0;
prior{im}.std=4;
prior{im}.norm=80;
prior{im}.m_true=2;
%prior{im}.seq_gibbs.step=1;
%prior{im}.seq_gibbs.i_update_step_max=2000;

% the intercept
im=im+1;
prior{im}.type='gaussian';
prior{im}.name='intercept';
prior{im}.m0=0;
prior{im}.std=15;
prior{im}.norm=10;
prior{im}.m_true=-10;
%prior{im}.seq_gibbs.step=1;
%prior{im}.seq_gibbs.i_update_step_max=2000;

prior=sippi_prior_init(prior);

%% Setup the forward model in the 'forward' structure
nd=40;
forward.x=sort(rand(1,nd)*20)';
forward.forward_function='sippi_forward_linefit';

%% Set up the 'data' structure
id=1;
m_true{1}=prior{1}.m_true;
m_true{2}=prior{2}.m_true;
d=sippi_forward_linefit(m_true,forward);
d_obs=d{1};
% Add noise top data
data{1}.d_std=10;
data{1}.d_obs=d_obs+randn(size(d_obs)).*data{1}.d_std;

%% Perform extended Metropolis sampling 
% set some MCMC options.
options.mcmc.nite=10000;
options.mcmc.i_sample=10;
options.mcmc.i_plot=1000;
options.txt='case_line_fit';
options=sippi_metropolis(data,prior,forward,options);


%% plot some stats
% get sample from posterior
m1_post=load(sprintf('%s%s%s_m1.asc',options.txt,filesep,options.txt))';
m1_post=m1_post(100:end);
m2_post=load(sprintf('%s%s%s_m2.asc',options.txt,filesep,options.txt))';
m2_post=m2_post(100:end);
% get sample prior
N=length(m1_post);
for i=1:1:N
    m=sippi_prior(prior);
    m1_prior(i)=m{1};
    m2_prior(i)=m{2};
end  

%% plot
try;close(1);end;figure(1);clf;set_paper('portrait');
x1=[0 20];
hold on
for i=ceil(linspace(1,N,100))
    plot(x1,m1_prior(i).*x1+m2_prior(i),'r-','linewidth',.8);
    plot(x1,m1_post(i).*x1+m2_post(i),'g-','linewidth',.8);
end
errorbar(forward.x,data{1}.d_obs,ones(size(forward.x)).*data{1}.d_std,'k.','linewidth',2)
hold off
box on
xlabel('X');ylabel('Y')
ppp(8,8,12)
print_mul('sippi_line_fit')

%%
try;close(2);end;figure(2);clf;set_paper('portrait');
plot(m1_prior,m2_prior,'r.','MarkerSize',6)
hold on
plot(m1_post,m2_post,'g.','MarkerSize',14)
plot(prior{1}.m_true,prior{2}.m_true,'k.','MarkerSize',24);

hold off
xlabel('Gradient')
ylabel('Intercept')
legend('A priori','A posteriori','true model')
ppp(8,8,12)
print_mul('sippi_line_fit_cross')

%%
for i=1:length(m1_post);
    m{1}=m1_post(i);
    m{2}=m2_post(i);
    d=sippi_forward_linefit(m,forward);
    logL(i)=sippi_likelihood(d,data);
end
maxL=max(logL),

%% REJECTION SAMPLING
clear options
options.mcmc.adaptive_rejection=1;
options.mcmc.nite=40000;
options=sippi_rejection(data,prior,forward,options);

% load the posterior sample
G_rej=load([options.txt,filesep,options.txt,'_m1.asc']);
I_rej=load([options.txt,filesep,options.txt,'_m2.asc']);

%
ii=1:min([length(m1_post) length(G_rej) 500]);
try;close(3);end;figure(3);clf;set_paper('portrait');
try
plot(m1_prior(ii),m2_prior(ii),'r.','MarkerSize',8)
hold on
plot(m1_post(ii),m2_post(ii),'g.','MarkerSize',6)
end
plot(G_rej(ii),I_rej(ii),'b.','MarkerSize',6)
try
plot(prior{1}.m_true,prior{2}.m_true,'k.','MarkerSize',24);
end
hold off
xlabel('Gradient')
ylabel('Intercept')
legend('\sigma','\rho_{Metropolis}','\rho_{Rejection}','true model')
ppp(8,8,12)
print_mul('sippi_line_fit_cross_rejection')

save sippi_line_fit

%% PLOT HIST OF logLikelihoods
hx=[-150:.5:-130];
h1=hist(logL,hx);h1=h1./sum(h1);
h2=hist(options.mcmc.logL,hx);h2=h2./sum(h2);
plot(hx,[h1;h2]','-');
b=bar(hx,[h1;h2]');
set(b(1),'EdgeColor','k');set(b(1),'FaceColor','g','LineWidth',.1,'BarWidth',1.4)
set(b(2),'EdgeColor','k');set(b(2),'FaceColor','b','LineWidth',.1,'BarWidth',1.4)
ppp(8,8,12,2,2)
legend('Metropolis','Rejection','Location','NorthWest')
xlabel('log(Likelihood)')
ylabel('Count (#)')
set(gca,'xlim',[hx(1) hx(length(hx))])
print_mul('sippi_line_fit_cross_metro_rejection')
%% LINEAR LEAST SQUARESSOLUTION
%sippi_least_squares()

