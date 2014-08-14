% prior_reals_gaussian Sampling two 1D GAUSSIAN type prior models

clear all;close all

ip=1;
prior{ip}.type='GAUSSIAN';
prior{ip}.name='Gaussian';
prior{ip}.m0=10;
prior{ip}.std=2;

ip=2;
prior{ip}.type='GAUSSIAN';
prior{ip}.name='Laplace';
prior{ip}.m0=10;
prior{ip}.std=2;
prior{ip}.norm=1;

ip=3;
prior{ip}.type='GAUSSIAN';
prior{ip}.name='Uniform';
prior{ip}.m0=10;
prior{ip}.std=2;
prior{ip}.norm=60;


ip=4;
prior{ip}.type='GAUSSIAN';
prior{ip}.name='TARGET';
d1=randn(100)+8;
d2=randn(200)+12;
prior{ip}.d_target=[d1(:);d2(:)];
prior{ip}.m0=0;


%%
sippi_plot_prior_sample(prior,1:length(prior),10000);
for ip=1:length(prior);
    figure(90+ip);print_mul(['prior_gaussian_1d_',prior{ip}.name])
end
return
%%
N=10000;
m_real1=ones(1,N);
m_real2=ones(1,N);
for i=1:N;
    [m,prior]=sippi_prior(prior);
    m_real1(i)=m{1};
    m_real2(i)=m{2};
end

%%
x=linspace(0,20,41);
h1=hist(m_real1,x);
h2=hist(m_real2,x);
figure(5);clf
plot(x,h1,'k-','linewidth',2)
hold on
plot(x,h2,'b-','linewidth',2)
hold off
xlabel('m')
ylabel('count #')
ppp(8,8,10)
legend('norm=60','norm=2')
print_mul('prior_reals_gaussian');

%% RANDOM WALK
prior{1}.seq_gibbs.step=0.1;
prior{2}.seq_gibbs.step=0.1;
N=10000;
m_real=ones(2,N);
prior=sippi_prior_init(prior);
[m,prior]=sippi_prior(prior);
for i=1:N;
    [m,prior]=sippi_prior(prior,m);
    m_real(1,i)=m{1};
    m_real(2,i)=m{2};
end

%%
np=300;
figure(6);clf
plot(m_real(:,1:np)','LineWidth',2)
legend('norm=60','norm=2')
xlabel('Iteration number')
ylabel('Realization value')
ppp(8,8,10)
print_mul('prior_reals_gaussian_seq_gibbs_progress');

%%
x=linspace(0,20,41);
h1=hist(m_real(1,:),x);
h2=hist(m_real(2,:),x);
figure(7);clf
plot(x,h1,'k-','linewidth',2)
hold on
plot(x,h2,'b-','linewidth',2)
hold off
xlabel('m')
ylabel('count #')
ppp(8,8,10)
legend('norm=60','norm=2')
print_mul('prior_reals_gaussian_seq_gibbs');

