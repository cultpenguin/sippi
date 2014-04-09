% prior_reals_gaussian Sampling two 1D GAUSSIAN type prior models

clear all;close all


im=1;
prior{im}.type='GAUSSIAN';
prior{im}.name='Gaussian';
prior{im}.m0=10;
prior{im}.std=2;

im=2;
prior{im}.type='GAUSSIAN';
prior{im}.name='Laplace';
prior{im}.m0=10;
prior{im}.std=2;
prior{im}.norm=1;

im=3;
prior{im}.type='GAUSSIAN';
prior{im}.name='Uniform';
prior{im}.m0=10;
prior{im}.std=2;
prior{im}.norm=60;

%%
sippi_plot_prior(prior,1:3,10000);
for im=1:length(prior);
    figure(90+im);print_mul(['prior_gaussian_1d_',prior{im}.name])
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

