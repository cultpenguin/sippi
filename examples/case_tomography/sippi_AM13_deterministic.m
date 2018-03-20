% sippi_AM13_least_squares 2D least squares inversion  
%
% Example of inverting 2D Arren�s tomographic data (AM13)
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
data{id}.d_std=0.*D.d_std+0.3;;
data{id}.Ct=D.Ct;
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

%% FORWARD
i_forward = 2;
forward.forward_function='sippi_forward_traveltime';
if i_forward==1;
    % STRAIGHT RAY FORWARD
    try; forward=rmfield(forward,'G');end
    forward.type='ray_2d';
    forward.linear=1;
    options.txt=[txt,'_',forward.type];
    
elseif i_forward==2;
    % LINEAR FAT FORWARD
    try; forward=rmfield(forward,'G');end
    forward.type='fat';
    forward.linear=1;
    options.txt=[txt,'_',forward.type];
    forward.freq=.1;
   
elseif i_forward==3;
    % LINEAR BORN FORWARD
    try; forward=rmfield(forward,'G');end
    forward.type='born';
    forward.linear=1;
    options.txt=[txt,'_',forward.type];
    forward.freq=.1;
end

%% LEAST SQUARES
[m_est,Cm_est,m_reals,options,data,prior,forward]=sippi_least_squares(data,prior,forward,options);


%% THIKONOV
for i=2;
    close all;
    options.lsq.tikhonov=i;
    options.lsq.interactive = 1;
    [m_est_tik,options_lsq]=sippi_tikhonov(data,prior,forward,options);   
end

%%
figure(10);subfigure(1,1,1);clf;;
d_norm_lsq=norm(options_lsq.lsq.d_obs-options_lsq.lsq.G*m_est{1}(:));

subplot(1,7,1)
imagesc(prior{1}.x,prior{1}.y,reshape(m_est{1},length(prior{1}.y),length(prior{1}.x)));
title('\bf{m}_{lsq}')
axis image;try;caxis(prior{1}.cax);end
%print_mul(sprintf('%s_mest_lsq',options.txt))
for i=1:6;
    subplot(1,7,i+1)
    m=m_reals(:,i+1);
    d_norm_real(i)=norm(options_lsq.lsq.d_obs-options_lsq.lsq.G*m);

    imagesc(prior{1}.x,prior{1}.y,reshape(m,length(prior{1}.y),length(prior{1}.x)));
    title(sprintf('\bf{m}^{*}_{%d}',i))

    axis image;try;caxis(prior{1}.cax);end
end
print_mul(sprintf('%s_mest_lsq_reals',options_lsq.txt))

%%
figure(3);
xl=xlim;
hold on;
for i=1:length(d_norm_real);
    plot(xl,[1 1].*d_norm_real(i),'r-','LineWidth',1);
end
plot(xl,[1 1].*d_norm_lsq,'g-','LineWidth',2);
hold off
axis([0 300 0 20])
print_mul(sprintf('%s_Lcurve_N_reals',options_lsq.txt))



