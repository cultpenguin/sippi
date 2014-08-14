clear all;%close all;

%% TARGET
N=10000;
prob_chan=0.25;
d1=randn(1,ceil(N*(1-prob_chan)))*.5+8.5;
d2=randn(1,ceil(N*(prob_chan)))*.5+11.5;
d_target=[d1(:);d2(:)];
[d_nscore,o_nscore]=nscore(d_target,1,1,min(d_target),max(d_target),0);
% UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION

%% PRIOR
ip=0;

ip=ip+1;
prior{ip}.type='fftma';
dx=.25;
prior{ip}.x=0:dx:10;
prior{ip}.y=0:dx:20;
prior{ip}.Cm='1 Sph(10,90,.25)';
prior{ip}.d_target=d_target;
m0_vert=linspace(0,2,length(prior{ip}.y))';
prior{ip}.m0=repmat(m0_vert,1,length(prior{ip}.x));

prior{ip}.cax=[8 14];

% CHOL
ip=ip+1;
prior{ip}=prior{ip-1};
prior{ip}.type='Cholesky';

%% VISIM
%ip=ip+1;
%prior{ip}=prior{ip-1};
%prior{ip}.type='visim';


[m,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m);
%% 
for ip=1:length(prior);
prior{ip}.seq_gibbs.step=0.01;
end
for i=1:30;
    [m,prior]=sippi_prior(prior,m);
    sippi_plot_prior(prior,m);
drawnow;
end
%%


return
sippi_plot_prior_sample(prior)