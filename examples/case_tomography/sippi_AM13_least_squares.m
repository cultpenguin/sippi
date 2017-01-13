% sippi_AM13_least_squares 2D least squares inversion  
%
% Example of inverting 2D Arrenæs tomographic data (AM13)
% using least squares inversion
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%


%%
%% Example of least squares tomography on AM13 dataset
%%
clear all;close all
D=load('AM13_data.mat');
txt='AM13';

forward.is_slowness=1; % use slowness parameterization

%% THE DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.i_use=10:10:702;

%% THE PRIOR
dx=0.2;im=1;

if forward.is_slowness==1
    % slowness
    prior{im}.type='FFTMA';
    prior{im}.name='Slowness (ns/m)';
    prior{im}.m0=7.0035;
    prior{im}.Va='0.7728 Exp(6)';
    prior{im}.x=[-1:dx:6];
    prior{im}.y=[0:dx:13];
    prior{im}.cax=1./[.18 .1];
else
    prior{im}.type='FFTMA';
    prior{im}.name='Velocity (m/ns)';
    prior{im}.m0=0.145;
    prior{im}.Va='.0003 Sph(6)';
    prior{im}.Va='.00001 Sph(6)';
    prior{im}.x=[-1:dx:6];
    prior{im}.y=[0:dx:13];
    prior{im}.cax=[-1 1].*.02+prior{im}.m0;
end
%% FORWARD MODEL
D=load('AM13_data.mat');
forward.sources=D.S;
forward.receivers=D.R;
%forward.is_slowness=1; % WE USE SLOWNESS PARAMETERIZAtiON

%% SOLVE LEAST SQUARES OPTIONS
options.lsq.type='lsq';
%lsq_type='visim';
%lsq_type='error_sim';

% set number of realization
options.lsq.n_reals=50;

% select a subset of data to consider
%data{1}.i_use=1:20:702;

%% STRAIGHT RAY FORWARD
forward.forward_function='sippi_forward_traveltime';
forward.type='ray';
forward.linear=1;
options.txt=[txt,'_',forward.type];

[m_est_1,Cm_est_1,m_reals_1,options_1,data_1,prior_1,forward_1]=sippi_least_squares(data,prior,forward,options);

%% LINEAR FAT FORWARD
try; forward=rmfield(forward,'G');end
forward.type='fat';
forward.linear=1;
options.txt=[txt,'_',forward.type];
forward.freq=.1;
[m_est_2,Cm_est_2,m_reals_2,options_2,data_2,prior_2,forward_2]=sippi_least_squares(data,prior,forward,options);

%% LINEAR BORN FORWARD
try; forward=rmfield(forward,'G');end
forward.type='born';
forward.linear=1;
options.txt=[txt,'_',forward.type];
forward.freq=.1;
[m_est_3,Cm_est_3,m_reals_3,options_3,data_3,prior_3,forward_3]=sippi_least_squares(data,prior,forward,options);

%% POST PLOTS
%sippi_plot_posterior_sample(options_1.txt)
%sippi_plot_posterior_sample(options_2.txt)
%sippi_plot_posterior_sample(options_3.txt)

%%
cax=fliplr(1./prior{1}.cax);
figure(3);clf;
subplot(1,3,1);
imagesc(prior{1}.x,prior{1}.y,1./m_est_1{1});caxis(cax);axis image
title('a) Ray kernel')
subplot(1,3,2);
imagesc(prior{1}.x,prior{1}.y,1./m_est_2{1});caxis(cax);axis image
title('a) Fat kernel')
subplot(1,3,3);
imagesc(prior{1}.x,prior{1}.y,1./m_est_3{1});caxis(cax);axis image
colorbar_shift;
title('a) Born kernel')
print_mul(sprintf('%s_compare_est',txt));

%%
prior=sippi_prior_init(prior);
figure(4);clf,set_paper('landscape')
nx=prior{1}.dim(1);
ny=prior{1}.dim(2);
for i=1:5;
    
    subplot(3,5,i);
    imagesc(prior{1}.x,prior{1}.y,reshape(1./m_reals_1(:,i),ny,nx))
    axis image;caxis(cax);
    
    subplot(3,5,i+5);
    imagesc(prior{1}.x,prior{1}.y,reshape(1./m_reals_2(:,i),ny,nx))
    axis image;caxis(cax);
    
    subplot(3,5,i+10);
    imagesc(prior{1}.x,prior{1}.y,reshape(1./m_reals_3(:,i),ny,nx))
    axis image;caxis(cax);
    
end
subplot(3,5,10);colorbar_shift;
try
    n_use=length(data{1}.i_use);
catch
    n_use=length(data{1}.d_obs);
end
print_mul(sprintf('%s_compare_reals_nd%d',txt,n_use));


