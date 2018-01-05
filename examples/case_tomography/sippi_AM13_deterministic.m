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

%% FORWARD
i_forward = 2;
forward.forward_function='sippi_forward_traveltime';
if i_forward==1;
    % STRAIGHT RAY FORWARD
    try; forward=rmfield(forward,'G');end
    forward.type='fat';
    forward.linear=1;
    options.txt=[txt,'_',forward.type];
    
elseif i_forward==2;
    % LINEAR FAT FORWARD
    try; forward=rmfield(forward,'G');end
    forward.type='fat';
    forward.linear=1;
    options.txt=[txt,'_',forward.type];
    forward.freq=.1;
    [m_est_2,Cm_est_2,m_reals_2,options_2,data_2,prior_2,forward_2]=sippi_least_squares(data,prior,forward,options);
    
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
sippi_thikonov(data,prior,forward)

