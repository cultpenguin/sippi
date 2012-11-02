% sippi_AM13_least_squares

%%
%% Example of least squares tomography on AM13 dataset
%%
clear all;close all
D=load('AM13_data.mat');
txt='AM13';

%% THE PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';
dx=0.2;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[0.1 0.18];
prior=sippi_prior_init(prior);
%% THE DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;
data{id}.Ct=0; % modelization error
data{id}.Ct=D.Ct; % modelization error in style of Cordua et al

%% FORWARD MODEL
D=load('AM13_data.mat');
forward.sources=D.S;
forward.receivers=D.R;

%% CONVERT MODEL TO SLOWNESS
forward.is_slowness=1; % WE USE SLOWNESS PARAMETERIZAtiON
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

forward.is_slowness=1; % WE USE SLOWNESS PARAMETERIZAtiON


%% SOLVE LEAST SQUARES OPTIONS
lsq_type='lsq';
%lsq_type='visim';
%lsq_type='error_sim';

% set number of realization
n_reals=15;

% select number of data to consider
%data{1}.i_use=1:50:702;

%% STRAIGHT RAY FORWARD
forward.forward_function='sippi_forward_traveltime';
forward.type='ray';
forward.linear=1;
options.txt=[txt,'_',forward.type];
[m_reals_1,m_est_1,Cm_est_1,options_1]=sippi_least_squares(data,prior,forward,n_reals,lsq_type,1,1,options);
%sippi_plot_posterior(options_1.txt)

%% LINEAR FAT FORWARD
forward.type='fat';
forward.linear=1;
options.txt=[txt,'_',forward.type];
forward.freq=.1;
[m_reals_2,m_est_2,Cm_est_2,options_2]=sippi_least_squares(data,prior,forward,n_reals,lsq_type,1,1,options);
sippi_plot_posterior(options_2.txt)
%% LINEAR BORN FORWARD
forward.type='born';
forward.linear=1;
options.txt=[txt,'_',forward.type];
forward.freq=.1;
[m_reals_3,m_est_3,Cm_est_3,options_3,forward]=sippi_least_squares(data,prior,forward,n_reals,lsq_type,1,1,options);
%sippi_plot_posterior(options_3.txt)

%%
cax=[.1 .18]
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
print_mul(sprintf('%s_compare_reals',txt));


