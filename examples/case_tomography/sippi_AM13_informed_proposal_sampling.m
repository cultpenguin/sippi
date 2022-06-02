% sippi_AM13_informed_proposal_sampling
%
% Example of inverting 2D Arrenæs tomographic data (AM13)
% using least squares inversion
%
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%


%%
%% Example of least squares tomography on AM13 dataset
%%
clear all;%close all
dx=0.4;
im=1;


D=load('AM13_data.mat');
use_sparse=1;
if use_sparse==1
    id_use = 1:4:size(D.S,1);
    D.S=D.S(id_use,:);
    D.R=D.R(id_use,:);
    D.d_obs=D.d_obs(id_use,:);
    D.d_std=D.d_std(id_use,:);%.*0+0.01;
    D.Ct=D.Ct(id_use,id_use);
end

txt='AM13';

%% FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.forward_function='sippi_forward_traveltime';
forward.is_slowness=1; % use slowness parameterization

%% THE PRIOR

if forward.is_slowness==1
    % slowness
    prior{im}.type='CHOLESKY';
    %prior{im}.type='FFTMA';
    prior{im}.name='Slowness (ns/m)';
    prior{im}.m0=7.0035;
    prior{im}.Cm='0.7728 Exp(6)';
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


%% THE DATA
useSynth=1;
if useSynth==1;
    rng(3)
    [m_ref,prior]=sippi_prior(prior);
    forward_ref = forward;
    forward_ref.type='eikonal';
    %forward_ref.type='ray';forward_ref.linear=1;
    %forward_ref.type='fat';forward_ref.linear=1;
    d_ref = sippi_forward(m_ref,forward_ref,prior);
    id=1;

    data{id}.d_ref=d_ref{1};
    data{id}.d_std=D.d_std;
    data{id}.d_noise=randn(size(data{id}.d_std)).*data{id}.d_std;
    data{id}.d_obs=data{id}.d_ref + data{id}.d_noise;    
else 
    id=1;
    data{id}.d_obs=D.d_obs;
    data{id}.d_std=D.d_std;
    %data{id}.i_use=10:10:702;
end


%% SOLVE LSQ - approximate
options.lsq.type='lsq';
options.lsq.compute_reals=0;
options.lsq.n_reals=50;
forward.forward_function='sippi_forward_traveltime';
forward.type='ray';forward.linear=1;
options.txt=[txt,'_',forward.type];

data_app=data;
data_app{1}.d_std = 4*data{1}.d_std;

[m_est,Cm_est,m_reals_1,options_1,data_1,prior,forward]=sippi_least_squares(data_app,prior,forward,options);
iCm=inv(Cm_est{1});
[~,~,LCm]=gaussian_simulation_cholesky(m_est{1},Cm_est{1},1);

%% McMC
clear m_post logL

usePert=2;
usePrior=1;

cax=[0.11, 0.17];
m_min = 4.5;
m_max = 9;
dm=m_max-m_min;
nm=prod(prior{1}.dim);

m_cur{1}=rand(size(m_ref{1}))*dm+m_min;
d_cur = sippi_forward(m_cur,forward_ref,prior);
L_cur = sippi_likelihood(d_cur,data);
if usePrior==1
    %% THIS IS ACTUALLY THE PROPOSAL DISTRIBUTION -- FIX
    dmm=m_cur{1}-m_est{1};
    rho_cur = -0.5*dmm(:)'*iCm*dmm(:);
else
    rho_cur = 1;
end    

k=0;
for i=1:50000

    % Perturb   
    if usePert==1;
        % Independent
        m_pro{1}=rand(size(m_ref{1}))*dm+m_min;
    elseif usePert==2;
        % Random perturbation
        i_pert=randomsample(nm,ceil(0.01*nm));
        m_pro{1}=m_cur{1};
        m_pro{1}(i_pert)=rand(1,length(i_pert))*dm+m_min;
    elseif usePert==3;
        % From app posterior
        %is_chol=0;
        %mm=gaussian_simulation_cholesky(m_est{1},Cm_est{1},1,is_chol);
        is_chol=1;
        mm=gaussian_simulation_cholesky(m_est{1},LCm,1,is_chol);
        m_pro{1}=reshape(mm,prior{1}.dim(2),prior{1}.dim(1));
    elseif usePert==4;
        % Informed Proposal
        m_tilde = m_est;
        d_tilde = sippi_forward(m_tilde,forward,prior);
        data_tilde = data;
        data_tilde{1}.d_obs=d_tilde{1};
        [m_app,Cm_app]=sippi_least_squares(data_tilde,prior,forward,options);

        m_std = std([m_app{1}(:),m_tilde{1}(:)]')';

        m_pro{1} = m_tilde{1} + reshape(randn(size(m_std)),prior{1}.dim(2),prior{1}.dim(1));
    end

    d_pro = sippi_forward(m_pro,forward_ref,prior);
    L_pro = sippi_likelihood(d_pro,data);
    if usePrior==1
        dmm=m_pro{1}-m_est{1};
        rho_pro = -0.5*dmm(:)'*iCm*dmm(:);
    else
        rho_pro=1;
    end
    
    % Accept
    P_acc =  exp( (L_pro+rho_pro)-(L_cur+rho_cur) );
    if P_acc > rand(1)
        d_cur = d_pro;
        m_cur = m_pro;
        L_cur = L_pro;
        rho_cur = rho_pro;
    end

    logL(i)=L_cur;
    logRho(i)=rho_cur;

    if mod(i,100)==0
        k=k+1;
        m_post(:,:,k)=m_cur{1};

        figure_focus(1);
        subplot(1,3,1);imagesc(prior{1}.x,prior{1}.y,1./m_ref{1})
        axis image;caxis(cax)
        subplot(1,3,2);imagesc(prior{1}.x,prior{1}.y,1./m_cur{1})
        axis image;caxis(cax)
        subplot(1,3,3);imagesc(prior{1}.x,prior{1}.y,1./etype(m_post))
        axis image;caxis(cax)
        drawnow;
        figure_focus(2);
        subplot(1,2,1);
        plot(1:i,logL(1:i),'k-');ylabel('log(L)')
        hold on
        plot(1:i,logRho(1:i),'r-');ylabel('log(rho)')
        hold off
        subplot(1,2,2);
        errorbar(data{1}.d_obs,2*data{1}.d_std,'k.');
        hold on
        plot(d_cur{1},'r*');
        hold off
        drawnow;
        
    end
end









return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if useSynth==1;
    subplot(1,4,1);
    imagesc(prior{1}.x,prior{1}.y,1./m_ref{1});caxis(cax);axis image
    title('Reference')
end
subplot(1,4,2);
imagesc(prior{1}.x,prior{1}.y,1./m_est_1{1});caxis(cax);axis image
title('a) Ray kernel')
subplot(1,4,3);
imagesc(prior{1}.x,prior{1}.y,1./m_est_2{1});caxis(cax);axis image
title('a) Fat kernel')
subplot(1,4,4);
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


