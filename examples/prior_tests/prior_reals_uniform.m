% prior_reals_uniform
prior{1}.type='uniform';
prior{1}.min=-1;
prior{1}.max=1;

sippi_plot_prior_sample(prior,1,1000)
print_mul('prior_uniform');


%% seq gibbs
prior{1}.seq_gibbs.step=0.2;
N=1000;
m_all=zeros(1,N);
[m,prior]=sippi_prior(prior);     
for i=1:1000;
     [m,prior]=sippi_prior(prior,m);
     m_all(i)=m{1};
     subplot(1,2,1);plot(1:i,m_all(1:i));
end
xlabel('Iteration #')
ylabel('m_1')

subplot(1,2,2);hist(m_all)
xlabel('m_1')
ylabel('frequency')
print_mul('prior_uniform_seqgibbs');


