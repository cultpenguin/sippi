% sippi_AM13_informed_proposal_sampling
%
% Comparing different approaches to sample the posterior distribution using
% * Linear Least Squqares (Tarantola and Valette, 1982)
% * The Metropolis-Hastings algortihm (MHA)
% * The exetended Metropolis algorithm (EMA)
% * The Metropolis algoriothm (MA)
% * Informed proposal algorithm (IPA)
% 
% 
% As an example travetime data from a GPR corss hole experiment at 
% Arrenæs, Denmark, is used 
% See http://dx.doi.org/10.1016/j.cageo.2012.10.001
%

%% Some general settinghs
clear all;
close all
useSynth=1;
dx=0.1;
d_std = 0.8;
%id_use = 1:1:702; % Selevt number of data to use
id_use = 1:1:702; % Use every 10th data

txt=sprintf('dx%g_std%g_S%d',dx*10,d_std*1000,useSynth);
cax=[0.11, 0.17];
cax_var = [0 0.7];

%% LOAD RAW DATA
D=load('AM13_data.mat');
D.d_std=D.d_std.*0+d_std;
if exist('id_use','var')
    %id_use = 1:4:size(D.S,1);
    D.S=D.S(id_use,:);
    D.R=D.R(id_use,:);
    D.d_obs=D.d_obs(id_use,:);
    D.d_std=D.d_std(id_use,:)
    D.Ct=D.Ct(id_use,id_use);
end
Nd=length(D.d_obs);

txt=sprintf('Nd%d_dx%g_std%g_S%d',Nd,dx*10,d_std*1000,useSynth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup PRIOR, DATA, and FORWARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FORWARD MODEL
% FULL FORWARD MODEL
forward.sources=D.S;
forward.receivers=D.R;
forward.type='eikonal';
%forward.type='ray';forward.linear=0; % bended ray
%forward.type='fat';forward.linear=0; % bended fat-ray
%forward.type='ray';forward.linear=1; % straight ray
%forward.type='fat';forward.linear=1; % straight fat-ray
forward.forward_function='sippi_forward_traveltime';
forward.is_slowness=1; % use slowness parameterization

% APPROXIMATE FORWARD
forward_app.sources=D.S;
forward_app.receivers=D.R;
forward_app.forward_function='sippi_forward_traveltime';
forward_app.type='ray';forward_app.linear=1; % straight ray
forward_app.type='fat';forward_app.linear=1; % straight fat
forward_app.is_slowness=1; % use slowness parameterization


%% THE PRIOR
im=1;
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
x=prior{1}.x;nx=length(x);
y=prior{1}.y;ny=length(y);

%% THE DATA
if useSynth==1;
    rng(3)
    [m_ref,prior]=sippi_prior(prior);
    d_ref = sippi_forward(m_ref,forward,prior);
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

%% 
m=sippi_prior(prior);
[d,forward,prior]=sippi_forward(m,forward,prior);
[d_app,forward_app]=sippi_forward(m,forward_app,prior);

if useSynth==1;
    d_ref = sippi_forward(m_ref,forward,prior);
    logL_ref = sippi_likelihood(d_ref,data);
    d_ref_app = sippi_forward(m_ref,forward_app,prior);
    logL_ref_app = sippi_likelihood(d_ref_app,data);
end

%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute modeling error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_model_err = 600;
[Ct,dt,dd]=sippi_compute_modelization_forward_error(forward,forward_app,prior,N_model_err);
data_me = data;
data_me{1}.Ct=Ct{1};
data_me{1}.dt=dt{1};
figure(11);clf
subplot(1,3,1);
histogram(dd{1});
xlabel('\Delta d')
ylabel('counts');title('Historgram of global modeling error')
subplot(1,3,2);
plot(data_me{1}.dt);grid on
xlabel('data #');ylabel('d_T')
subplot(1,3,3);imagesc(data_me{1}.Ct);axis image;colorbar
xlabel('data #');ylabel('data #');title('C_T')
sgtitle(sprintf('Estimated \\theta(d|m), N=%d',N_model_err))
print_mul(sprintf('%s_modelerror',txt))

%% SOLVE LSQ - approximate
options.lsq.type='lsq';
options.lsq.compute_reals=1;
options.lsq.n_reals=150;

[m_lsq,Cm_lsq,m_reals_lsq,o_lsq]=sippi_least_squares(data,prior,forward_app,options);
[m_lsq_me,Cm_lsq_me,m_reals_lsq_me,o_lsq_me]=sippi_least_squares(data_me,prior,forward_app,options);

% precalc inv(Cm) and chol(Cm)
iCm_lsq{1}=inv(Cm_lsq{1});
iCm_lsq_me{1}=inv(Cm_lsq_me{1});
[~,~,L_lsq]=gaussian_simulation_cholesky(m_lsq{1},Cm_lsq{1});
[~,~,L_lsq_me]=gaussian_simulation_cholesky(m_lsq_me{1},Cm_lsq_me{1});

%% compute log(L)
for i=1:size(m_reals_lsq,2)
    progress_txt(i,size(m_reals_lsq,2),'logL computation')
    m_test{1} = reshape(m_reals_lsq(:,i),ny,nx);
    d_test = sippi_forward(m_test,forward,prior);
    logL_lsq(i) = sippi_likelihood(d_test,data);

    m_test{1} = reshape(m_reals_lsq_me(:,i),ny,nx);
    d_test = sippi_forward(m_test,forward,prior);
    logL_lsq_me(i) = sippi_likelihood(d_test,data);

end

%% 
% PLOT LSQ RESULTS
figure(12);clf;
if useSynth==1;
    subplot(2,3,1);
    imagesc(x,y,1./m_ref{1});axis image;caxis(cax);
    title('m_{ref}')
end

subplot(2,3,2);
imagesc(x,y,1./m_lsq{1});axis image;caxis(cax);
title('m_{lsq}')
subplot(2,3,5);
imagesc(x,y,reshape(diag(Cm_lsq{1}),ny,nx));axis image;caxis(cax_var)
title('Cm_{lsq}')

subplot(2,3,3);
imagesc(x,y,1./m_lsq_me{1});axis image;caxis(cax)
title('m_{lsq_{me}}')
subplot(2,3,6);
imagesc(x,y,reshape(diag(Cm_lsq_me{1}),ny,nx));axis image;caxis(cax_var)
title('Cm_{lsq_{me}}')
sgtitle('Least Sqaures')
print_mul(sprintf('%s_lsq_est',txt))

% plot posterior realizations!
%print_mul(sprintf('%s_lsq_real',txt))

% PLOT likelihood
figure(13);clf
plot(logL_lsq,'r--')
hold on
plot(logL_lsq_me,'b:')
plot(xlim,[logL_ref logL_ref],'k-','LineWidth',2)
plot(xlim,[logL_ref_app logL_ref_app],'k-.','LineWidth',2)
plot(xlim,[-1 -1].*Nd/2,'k--','LineWidth',2)
hold off
grid on
legend('Lsq','Lsq_{me}','logL(g(m_{ref}))','logL(g_{app}(m_{ref}))','-N/2')
print_mul(sprintf('%s_lsq_logl',txt))

%%%



%% McMC
prior{1}.seq_gibbs.step=.5;

iCm_prior = inv(prior{1}.Cmat);
L_prior=chol(prior{1}.Cmat,'lower');
m0_prior = prior{1}.m0;
clear  m_post

usePert=3; usePrior=4; % Informed proposal
%usePert=3; usePrior=3; % Extended Metropolis
useME=1;

m_min = 4.5;
m_max = 9;
dm=m_max-m_min;
nm=prod(prior{1}.dim);

% start model
if usePrior==1
    m_cur{1}=rand(size(m_ref{1}))*dm+m_min;
elseif usePrior==2
    m_cur = sippi_prior(prior);
elseif usePrior==3
    m_cur{1}=reshape(gaussian_simulation_cholesky(m0_prior,L_prior,1,1),ny,nx);
elseif usePrior==4
    m_cur{1}=reshape(gaussian_simulation_cholesky(m_lsq_me{1},L_lsq_me,1,1),ny,nx);
end

if useME==0;
    dmm=m_cur{1}-m_lsq{1};
    q_cur = -0.5*dmm(:)'*iCm_lsq{1}*dmm(:);
else
    dmm=m_cur{1}-m_lsq_me{1};
    q_cur = -0.5*dmm(:)'*iCm_lsq_me{1}*dmm(:);
end

dmm=m_cur{1}-m0_prior;
rho_cur = -0.5*dmm(:)'*iCm_prior*dmm(:);

d_cur = sippi_forward(m_cur,forward,prior);
L_cur = sippi_likelihood(d_cur,data);
%%

txt_sampling=sprintf('%s_PE%d_PR%d',txt,usePert,usePrior);
k=0;
for i=1:20000

    if i<1000, 
        step=0.5;
        prior{1}.seq_gibbs.step=step;
    elseif i<2000;
        step=0.2;
        prior{1}.seq_gibbs.step=step;
    else
        step=0.1;
        prior{1}.seq_gibbs.step=step;
    end
    % Perturb   
    if usePert==1;
        % Independent
        m_pro{1}=rand(size(m_ref{1}))*dm+m_min;
    elseif usePert==2;
        % Random perturbation
        i_pert=randomsample(nm,ceil(0.01*nm));
        m_pro{1}=m_cur{1};
        m_pro{1}(i_pert)=rand(1,length(i_pert))*dm+m_min;
    elseif usePert==3; % Random walk

        if usePrior==1
            m_cur{1}=rand(size(m_ref{1}))*dm+m_min;
        elseif usePrior==2
            m_pro = sippi_prior(prior,m_cur);
        elseif usePrior==3
            % Sequential Gibbs
            m=reshape(gaussian_simulation_cholesky(m0_prior,L_prior,1,1),ny,nx);
            step = prior{1}.seq_gibbs.step;
            m_pro{1} = m0_prior+cos(step*pi/2)*(m_cur{1}-m0_prior) + sin(step*pi/2)*(m-m0_prior);
        elseif usePrior==4
            % Use Proposal from linear inversion!!!
            m=reshape(gaussian_simulation_cholesky(m_lsq_me{1},L_lsq_me,1,1),ny,nx);
            step = prior{1}.seq_gibbs.step;
            m_pro{1} = m_lsq_me{1} + cos(step*pi/2)*(m_cur{1}-m_lsq_me{1}) + sin(step*pi/2)*(m-m_lsq_me{1});
        end

    elseif usePert==4;
        % Informed Proposal
        m_tilde = m_est;
        d_tilde = sippi_forward(m_tilde,forward,prior);
        data_tilde = data;
        data_tilde{1}.d_obs=d_tilde{1};
        [m_app,Cm_app]=sippi_least_squares(data_tilde,prior,forward_app,options);

        m_std = std([m_app{1}(:),m_tilde{1}(:)]')';

        m_pro{1} = m_tilde{1} + reshape(randn(size(m_std)),prior{1}.dim(2),prior{1}.dim(1));
    end

    dmm=m_pro{1}-m0_prior;
    rho_pro = -0.5*dmm(:)'*iCm_prior*dmm(:);

    if useME==0;
        dmm=m_pro{1}-m_lsq{1};
        q_pro = -0.5*dmm(:)'*iCm_lsq{1}*dmm(:);
    else
        dmm=m_pro{1}-m_lsq_me{1};
        q_pro = -0.5*dmm(:)'*iCm_lsq_me{1}*dmm(:);
    end

    d_pro = sippi_forward(m_pro,forward,prior);
    L_pro = sippi_likelihood(d_pro,data);
    

    % Accept
    if usePrior==3
        % extended type
        P_acc =  exp( (L_pro)-(L_cur) );
    elseif usePrior ==1, 
        P_acc =  exp( (L_pro+rho_pro)-(L_cur+rho_cur) );
    else
        P_acc =  exp( (L_pro+rho_pro+q_cur)-(L_cur+rho_cur+q_pro) );
    end
    if P_acc > rand(1)
        d_cur = d_pro;
        m_cur = m_pro;
        L_cur = L_pro;
        rho_cur = rho_pro;
        q_cur = q_pro;
    end

    logL(i)=L_cur;
    logRho(i)=rho_cur;
    logq(i)=q_cur;

    if mod(i,100)==0
        disp(sprintf('i=%d, step=%g',i,step));
    
        k=k+1;
        m_post(:,:,k)=m_cur{1};

        figure_focus(1);
        subplot(1,4,1);imagesc(prior{1}.x,prior{1}.y,1./m_ref{1})
        axis image;caxis(cax)
        subplot(1,4,2);imagesc(prior{1}.x,prior{1}.y,1./m_pro{1})
        axis image;caxis(cax);title('m_{pro}')
        subplot(1,4,3);imagesc(prior{1}.x,prior{1}.y,1./m_cur{1})
        axis image;caxis(cax);title('m_{cur}')
        subplot(1,4,4);imagesc(prior{1}.x,prior{1}.y,1./etype(m_post))
        axis image;caxis(cax)
        drawnow;
        figure_focus(2);
        subplot(1,2,1);
        plot(1:i,logL(1:i),'k-');ylabel('log(L(m))')
        hold on
        plot(1:i,logRho(1:i),'r-');ylabel('logL(rho(m))')
        plot(1:i,logq(1:i),'b-');ylabel('logL(q(m))')
        plot(xlim,[logL_ref logL_ref],'k--','LineWidth',2)
        plot(xlim,[-1 -1].*Nd/2,'k:','LineWidth',2)
        plot(xlim,[-1 -1].*(nx*ny)/2,'r:','LineWidth',2)
        hold off
        legend('log(L(m))','log(rho(m))','log(q(m))','log(L(m_ref))','-Nd/2','-Nm/2','Location','SouthEast')
        grid on
        subplot(1,2,2);
        errorbar(data{1}.d_obs,2*data{1}.d_std,'k.');
        hold on
        plot(d_cur{1},'r*');
        hold off
        drawnow;
        
    end
end



%% PLOT RESULTS
% PLOT LSQ RESULTS
figure(12);clf;
if useSynth==1;
    subplot(2,4,1);
    imagesc(x,y,1./m_ref{1});axis image;caxis(cax);
    title('m_{ref}')
end

subplot(2,4,2);
imagesc(x,y,1./m_lsq{1});axis image;caxis(cax);
title('m_{lsq}')
subplot(2,4,6);
imagesc(x,y,reshape(diag(Cm_lsq{1}),ny,nx));axis image;caxis(cax_var)
title('Cm_{lsq}')

subplot(2,4,3);
imagesc(x,y,1./m_lsq_me{1});axis image;caxis(cax)
title('m_{lsq_{me}}')
subplot(2,4,7);
imagesc(x,y,reshape(diag(Cm_lsq_me{1}),ny,nx));axis image;caxis(cax_var)
title('Cm_{lsq_{me}}')

% SIM
[em,ev]=etype(m_post(:,:,10:end));
subplot(2,4,4);
imagesc(x,y,1./em);axis image;caxis(cax)
title('m_{sim}')
subplot(2,4,8);
imagesc(x,y,ev);axis image;caxis(cax_var)
title('Cm_{sim}')

sgtitle('Posterior mean and std')
print_mul(sprintf('%s_est_compare',txt_sampling))


%%
figure(14);clf
h=histogram(logL_lsq);
dx=diff(h.BinEdges);dx=dx(1);
h_lsq_Values = h.Values./(dx*sum(h.Values));
h_lsq_x = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;

h=histogram(logL_lsq_me);
dx=diff(h.BinEdges);dx=dx(1);
h_lsq_me_Values = h.Values./(dx*sum(h.Values));
h_lsq_me_x = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;

h=histogram(logL);
dx=diff(h.BinEdges);dx=dx(1);
h_Values = h.Values./(dx*sum(h.Values));
h_x = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;

clf
plot(h_x,h_Values,'k-');
hold on
plot(h_lsq_me_x,h_lsq_me_Values);
plot(h_lsq_x,h_lsq_Values,'g-');
hold off
legend('LSQ','LSQ_{ME}','SAMPLING')
xlabel('Log(L)')

print_mul(sprintf('%s_logL_compare',txt_sampling))


save(txt_sampling)

% plot posterior realizations!
%print_mul(sprintf('%s_lsq_real',txt))





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


