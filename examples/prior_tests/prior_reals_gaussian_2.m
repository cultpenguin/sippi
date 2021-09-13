% prior_reals_gaussian Sampling two 1D GAUSSIAN type prior models

clear all;close all

ip=1;
prior{ip}.type='GAUSSIAN';
prior{ip}.name='Gaussian - norm 60';
prior{ip}.m0=5;
prior{ip}.std=2;
prior{ip}.norm=60;

ip=2;
prior{ip}.type='GAUSSIAN';
prior{ip}.name='Gaussian - target';
d1=randn(100)+8;
d2=randn(200)+12;
prior{ip}.d_target=[d1(:);d2(:)];
prior{ip}.m0=0;

%%
N=5000;
m_real1=ones(1,N);
m_real2=ones(1,N);
for i=1:N;
    [m,prior]=sippi_prior(prior);
    m_real1(i)=m{1};
    m_real2(i)=m{2};
end

%
x=linspace(0,20,41);
h1=hist(m_real1,x);
h2=hist(m_real2,x);
figure(5);clf
plot(x,h1,'r-','linewidth',2)
hold on
plot(x,h2,'b-','linewidth',2)
hold off
xlabel('m')
ylabel('count #')
%ppp(8,8,10)
legend(prior{1}.name,prior{2}.name)
sippi_plot_set_axis;
print_mul('prior_reals_gaussian_2');

%% RANDOM WALK
prior{1}.seq_gibbs.step=0.1;
prior{2}.seq_gibbs.step=0.1;
m_real=ones(2,N);
prior=sippi_prior_init(prior);
[m,prior]=sippi_prior(prior);
for i=1:N;
    [m,prior]=sippi_prior(prior,m);
    m_real(1,i)=m{1};
    m_real(2,i)=m{2};
end
%%
figure
subplot(2,1,1);
plot(m_real(:,:)','.')
xlabel('Iteration number')
ylabel('Vaule')
legend(prior{1}.name,prior{2}.name)
sippi_plot_set_axis;
subplot(2,2,3);
hist(m_real(1,:),40);
ylabel('pdf, prior 1')
sippi_plot_set_axis;
subplot(2,2,4);
hist(m_real(2,:),40);
ylabel('pdf, prior 2')
sippi_plot_set_axis;
disp('A')
print_mul('prior_reals_gaussian_2_seq_gibbs');
