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


%% THE DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
%data{id}.Ct=0; % modelization error
%data{id}.Ct=D.Ct; % modelization error in style of Cordua et al (2008; 2009)

%% THE PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
prior{im}.Va='.00001 Sph(6)';
dx=0.2;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[0.12 0.16];

%% FORWARD MODEL
D=load('AM13_data.mat');
forward.sources=D.S;
forward.receivers=D.R;
forward.is_slowness=1; % WE USE SLOWNESS PARAMETERIZAtiON

%% CONVERT PRIOR MODEL TO SLOWNESS
if forward.is_slowness==1
    if isstr(prior{im}.Va);
        prior{im}.Va=deformat_variogram(prior{im}.Va);
    end
    Va=prior{im}.Va;
    v_var=sum([Va.par1]);
    v_sim=randn(1,10000).*sqrt(v_var);
    slo_var=var(1./(v_sim+mean(prior{im}.m0(:))));
    slo_mean=mean(1./(v_sim+mean(prior{im}.m0(:))));
    Va_slo=slo_var*[Va.par1]./v_var;
    for i=1:length(Va);
        Va(i).par1=Va_slo(i);
    end
    prior{im}.Va=Va;
    prior{im}.m0=slo_mean;
    prior{1}.cax=fliplr(1./prior{1}.cax);
end

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

[options_1,data_1,prior_1,forward_1,m_reals_1,m_est_1,Cm_est_1]=sippi_least_squares(data,prior,forward,options);
sippi_plot_posterior_sample(options_1.txt)
%% LINEAR FAT FORWARD
try; forward=rmfield(forward,'G');end
forward.type='fat';
forward.linear=1;
options.txt=[txt,'_',forward.type];
forward.freq=.1;
[options_2,data_2,prior_2,forward_2,m_reals_2,m_est_2,Cm_est_2]=sippi_least_squares(data,prior,forward,options);
sippi_plot_posterior_sample(options_2.txt)

%% LINEAR BORN FORWARD
try; forward=rmfield(forward,'G');end
forward.type='born';
forward.linear=1;
options.txt=[txt,'_',forward.type];
forward.freq=.1;
[options_3,data_3,prior_3,forward_3,m_reals_3,m_est_3,Cm_est_3]=sippi_least_squares(data,prior,forward,options);
sippi_plot_posterior_sample(options_3.txt)

%%
cax=fliplr(1./prior{1}.cax);
figure(3);clf;
subplot(1,3,1);
imagesc(prior{1}.x,prior{1}.y,1./m_est_1);caxis(cax);axis image
title('a) Ray kernel')
subplot(1,3,2);
imagesc(prior{1}.x,prior{1}.y,1./m_est_2);caxis(cax);axis image
title('a) Fat kernel')
subplot(1,3,3);
imagesc(prior{1}.x,prior{1}.y,1./m_est_3);caxis(cax);axis image
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
print_mul(sprintf('%s_compare_reals_nd%d',txt,length(data{1}.i_use)));


