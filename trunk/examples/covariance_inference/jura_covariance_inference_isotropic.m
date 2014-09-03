% jura_covariance_inference
%
% Example of inferring properties of a Gaussian model from point data
%

%% LOAD THE JURA DATA
clear all;close all
[d_prediction,d_transect,d_validation,h_prediction,h_transect,h_validation,x,y,pos_est]=jura;
ix=1;
iy=2;
id=6;

% get the position of the data
pos_known=[d_prediction(:,[ix iy])];  

% perform normal score transformation of tha original data
[d,o_nscore]=nscore(d_prediction(:,id));
h_tit=h_prediction{id};

%% SETUP A PRIORI MODEL / ONLY RANGE --> ISTROPIC COVARIANCE MODEL
im=0;
% A close to uniform distribution of the range, U[0;3].
im=im+1;
prior{im}.type='gaussian';
prior{im}.name='range_1';
prior{im}.min=0.01;
prior{im}.max=3;
prior{im}.norm=100;


%% DATA
data{1}.d_obs=d; % observed data
data{1}.d_std=0;.5; % uncertainty of observed data (in form of standard deviation of the noise)
%data{1}.i_use=1:1:30;

%% FORWARD
forward.forward_function='sippi_forward_covariance_inference';
forward.point_support=1;
forward.pos_known=pos_known;
forward.stabilize=0;
% initial choice of N(m0,Cm), mean and sill are 0, and 1, due
% due to normal score
forward.m0=0;
forward.Cm='1 Sph(2)';


%% METROPOLIS SAMPLING
options.mcmc.nite=5000;
options.mcmc.i_plot=100;
options.mcmc.i_sample=10;
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)

sippi_plot_posterior(options.txt);

return

%% CLASSICAL 'EXPERIMENTAL' COVARIANCE INFERENCE
figure(1);clf;
[gamma_exp1,h_exp1]=semivar_exp(pos_known,d,20);
[gamma_exp_ang,h_exp,h_anh]=semivar_exp(pos_known,d,20,4);
plot(h_exp1,gamma_exp1,'k-','LineWidth',5);
hold on
plot(h_exp,gamma_exp_ang,'-');
hold off



%% TEST Forward model
m=sippi_prior(prior);
m{1}=1;
[d,forward,prior,data]=sippi_forward(m,forward,prior,data);
[logL,L,data]=sippi_likelihood(d,data);

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
title(h_prediction(idata))

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

    