% jura_covariance_inference
%
% Example of inferring properties of a Gaussian model from point data
%

%% LOAD THE JURA DATA
clear all;close all
[d_prediction,d_transect,d_validation,h_prediction,h_transect,h_validation,x,y,pos_est]=jura;
ix=1;
iy=2;
id=6;         % select the id of the property from the Jura data to use
use_nscore=0; % Use normal score=

name=h_prediction{id};
% get the position of the data
pos_known=[d_prediction(:,[ix iy])];  

if use_nscore==0;
    d_obs=d_prediction(:,id);
else
    % perform normal score transformation of tha original data
    [d_obs,o_nscore]=nscore(d_prediction(:,id));
    h_tit=h_prediction{id};
    name=[name,'_nscore'];
end
%% SETUP A PRIORI MODEL / ONLY RANGE --> ISTROPIC COVARIANCE MODEL
im=0;
% A close to uniform distribution of the range, U[0;3].
im=im+1;
prior{im}.type='uniform';
prior{im}.name='range_1';
prior{im}.min=0;
prior{im}.max=3;

im=im+1;
prior{im}.type='uniform';
prior{im}.name='range_2';
prior{im}.min=0;
prior{im}.max=3;

im=im+1;
prior{im}.type='uniform';
prior{im}.name='ang_1';
prior{im}.min=0;
prior{im}.max=90;

im=im+1;
prior{im}.type='uniform';
prior{im}.name='nugget_fraction';
prior{im}.min=0;
prior{im}.max=1;


%% DATA
data{1}.d_obs=d_obs; % observed data
data{1}.d_std=0.1*std(d_obs);.4; % uncertainty of observed data (in form of standard deviation of the noise)
%data{1}.i_use=1:1:30;

%% FORWARD
forward.forward_function='sippi_forward_covariance_inference';
forward.point_support=1;
forward.pos_known=pos_known;
forward.stabilize=0;
% initial choice of N(m0,Cm), mean and sill are 0, and 1, due
% due to normal score
forward.m0=mean(d_obs);
forward.Cm=sprintf('%3.1f Sph(2)',var(d));

%%
[m,prior]=sippi_prior(prior);
[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
[L,tmp,data]=sippi_likelihood(d,data);

%% METROPOLIS SAMPLING
for ip=1:length(prior);
    prior{ip}.seq_gibbs.i_update_step_max=3000;
end
options.plot.hardcopy_types=0; % no hardcopy
options.mcmc.nite=10000;
options.mcmc.i_plot=1000;
options.mcmc.i_sample=25;
options.txt=name;
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);

sippi_plot_prior(options.txt);
sippi_plot_posterior(options.txt);

return

%% CLASSICAL 'EXPERIMENTAL' COVARIANCE INFERENCE
figure(1);clf;
[gamma_exp1,h_exp1]=semivar_exp(pos_known,d_obs,20);
[gamma_exp_ang,h_exp,h_anh]=semivar_exp(pos_known,d_obs,20,4);
plot(h_exp1,gamma_exp1,'k-','LineWidth',5);
hold on
plot(h_exp,gamma_exp_ang,'-');
hold off



%% TEST Forward model
m=sippi_prior(prior);
m{1}=2;
[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
[logL,L,data]=sippi_likelihood(d,data)

%%
r1=linspace(prior{1}.min,prior{1}.max,51);
%r1=linspace(0,5,41);
for i=1:length(r1);
    m{1}=r1(i);
    [d,forward,prior,data]=sippi_forward(m,forward,prior,data);
    [ll(i),L,data]=sippi_likelihood(d,data);
end
%ll=ll-max(ll);
figure(2);clf;
plot(r1,exp(ll));xlabel('range_1');ylabel('logL(range_1)');
title(h_prediction(id))

%
iop=find(ll==max(ll));iop=iop(1);
m{1}=r1(iop);
[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
[logL,L,data]=sippi_likelihood(d,data);

%%
[gamma_exp,h]=semivar_exp(forward.pos_known,data{1}.d_obs,[0:.2:4]);
gamma_synth=semivar_synth(forward.Va,h);
figure(31);clf;
plot(h,gamma_exp,'k-',h,gamma_synth,'r*')

    
