% sippi_AM13_nonGaussianModelingError
%
%
% Ex:
% clear all;
% use_prior=3;
% use_metropolis=1;
% use_rejection=0;
%
% clear all; use_prior=3; use_metropolis=1;use_rejection=0;n_ite=1000;sippi_AM13
%
% clear all;close all;rseed=1;n_ite=50000;use_prior=2;doAnneal=1;Nme=2000;sippi_AM13_nonGaussianModerlingError

if exist('TEST','var')
    if TEST
        close all,clear all;
        use_forward=4;
        use_forward_ref=5;
        use_prior=1;
        dx=0.2;
        rseed=1;
        useZOP=0;
        i_use=10:10:702;
        n_ite=50000;
        sippi_AM13_nonGaussianModerlingError
        %use_prior=2;
        %use_forward=1;
        %use_forward_ref=6;
        %rseed=1;
        %n_ite=500000;
        %doAnneal=1;
        %Nme=2000;
    end
end

% Load data from Arrenæs
D=load('AM13_data.mat');

%% Make some choices
if ~exist('use_prior','var')
    use_prior=1; % Gaussian
    % use_prior=2; % Gaussian with bimodal target distribution
    % use_prior=3; % Gaussian with uniform target distribution
    % use_prior=4; % Pnlurigaussian
    % use_prior=5; % Gaussian with variable covariance parameters
    % use_prior=6; % Matern type covariance with varying nu parameter
end
if ~exist('use_forward','var')
    %use_forward=1; % ray_2d - linear straight ray
    % use_forward=2; % linear ray (using eikonal)
    % use_forward=3; % bended ray
    use_forward=4; % linear fat
    % use_forward=5; % bended fat
    % use_forward=6; % eikonal
    % use_forward=7; % waveform FD + first arriavle
end

if ~exist('use_forward_ref','var')
    % use_forward_ref=1; % ray_2d - linear straight ray
    % use_forward_ref=2; % linear ray (using eikonal)
    % use_forward_ref=3; % bended ray
    % use_forward_ref=4; % linear fat
    use_forward_ref=5; % bended fat
    % use_forward_ref=6; % eikonal
    % use_forward_ref=7; % waveform FD + first arriavle
end


if ~exist('use_metropolis','var')
    use_metropolis=1;
end
if ~exist('use_rejection','var')
    use_rejection=0;
end

if ~exist('use_reference','var')
    use_reference=1;
end

if ~exist('use_correlated_noise','var')
    use_correlated_noise=0;
end
if ~exist('n_ite','var')
    n_ite=200000;
end
if ~exist('rseed','var')
    rseed=round(rand*1000000);
end
if ~exist('n_reals_out','var')
    n_reals_out=200;
end

if ~exist('doAnneal','var')
    doAnneal=0;
end
if ~exist('doTempering','var')
    doTempering=0;
end
if ~exist('dx','var')
    dx=0.2;
end
if ~exist('cov_range','var')
    cov_range=6;
end
if ~exist('Nme','var')
    Nme=1000;
end

if ~exist('i_use','var')
    if ~exist('useZOP','var'), useZOP=1;end
    if useZOP==1
        z_pos=unique(D.S(:,2));
        i_use=[];
        for iz=1:length(z_pos);
            for i=1:length(D.S(:,2));
                if (D.S(i,2)==z_pos(iz))&(D.R(i,2)==z_pos(iz))
                    i_use=[i_use,i];
                end
            end
        end
    else
        i_use=1:1:702;
        i_use=find(D.S(:,2)==6);
    end
end

if exist('plot_posterior','var');
    plot_posterior_sample=plot_posterior;
    plot_posterior_data=plot_posterior;
    plot_posterior_loglikelihood=plot_posterior;
    plot_posterior_2d_marg=plot_posterior;
end

if ~exist('plot_posterior_sample','var');    plot_posterior_sample=1;end
if ~exist('plot_posterior_data','var');    plot_posterior_data=1;end
if ~exist('plot_posterior_loglikelihood','var');    plot_posterior_loglikelihood=1;end
if ~exist('plot_posterior_2d_marg','var');    plot_posterior_2d_marg=1;end

%% SETUP DATA, PRIOR and FORWARD

options.txt=sprintf('AM13_f%d_f%d_P%d_nd%d_nit%d',use_forward,use_forward_ref,use_prior,length(i_use),n_ite);

%% SETUP DIFFERENT PRIOR STRUCTURES
% define some standard values
m0=0.145;
Va=sprintf('.0003 Sph(%g,90,.3)',cov_range);
%Va=sprintf('.0003 Sph(%g,45,.3)',cov_range);

% some parameters needed by all a priori types
np=2;
prior_ref{1}.name='Velocity (m/ns)';
prior_ref{1}.x=[(0-np*dx):dx:(5+np*dx)];
prior_ref{1}.y=[(1-np*dx):dx:(12+np*dx)];
prior_ref{1}.cax=[.10 .17];

% define a number of a priori models
im_all=0;

%% GAUSSIAN
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian';
prior_all{im_all}{im}.type='FFTMA';
prior_all{im_all}{im}.m0=m0;
prior_all{im_all}{im}.Va=Va;

%% GAUSSIAN WITH BIMODAL TARGET DISTRIBUTION
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian Bimodal';
prior_all{im_all}{im}.type='FFTMA';
prior_all{im_all}{im}.Cm=Va;
% bimodal distribution
N=10000;
prob_chan=0.5;
dd=.015;
d1=randn(1,ceil(N*(1-prob_chan)))*.0025+0.145-dd;  %0.1125;
d2=randn(1,ceil(N*(prob_chan)))*.0025+0.145+dd; %0.155;
d_target=[d1(:);d2(:)];
prior_all{im_all}{im}.d_target=d_target;

%% GAUSSIAN WITH UNIFORM TARGET DISTRIBUTION
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian Uniform';
prior_all{im_all}{im}.type='FFTMA';
prior_all{im_all}{im}.Cm=Va;
N=10000;
prior_all{im_all}{im}.d_target=(rand(1,N)-.5)*.09+m0;

%% PLURIGAUSSIAN
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Plurigaussian';
prior_all{im_all}{im}.type='plurigaussian';
% 1D plurigaussian (truncated)
prior_all{im_all}{im}.pg_prior{1}.Cm=' 1 Gau(5,90,.5)';
prior_all{im_all}{im}.pg_map=[0.11 .11 .13 .13 .15 .15 .17 .13 .17];
% 2D plurigaussian
prior_all{im_all}{im}.pg_prior{2}.Cm=' 1 Sph(10,45,.2)';
prior_all{im_all}{im}.pg_map=[0.11 .11 .17 ; .13 .15 .17 ; .17 .13 .17];

%% GAUSSIAN with variable covariance model parameters
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian - variable pars';

prior_all{im_all}{im}.type='fftma';
prior_all{im_all}{im}.Va=Va;
prior_all{im_all}{im}.m0=m0;
i_master=im;

% range - horizontal
im=im+1;
prior_all{im_all}{im}.type='uniform';
prior_all{im_all}{im}.name='range_1';
prior_all{im_all}{im}.min=1;
prior_all{im_all}{im}.max=6;
prior_all{im_all}{im}.prior_master=i_master;

% range - vertical
im=im+1;
prior_all{im_all}{im}=prior_all{im_all}{im-1};
prior_all{im_all}{im}.name='range_2';

% rotation
im=im+1;
prior_all{im_all}{im}.type='gaussian';
prior_all{im_all}{im}.name='ang_1';
prior_all{im_all}{im}.m0=90;
prior_all{im_all}{im}.std=10;
prior_all{im_all}{im}.norm=2;
prior_all{im_all}{im}.prior_master=i_master;

% m0
im=im+1;
prior_all{im_all}{im}.type='uniform';
prior_all{im_all}{im}.name='m0';
prior_all{im_all}{im}.min=0.10;
prior_all{im_all}{im}.max=0.17;
prior_all{im_all}{im}.prior_master=i_master;

% sill
im=im+1;
prior_all{im_all}{im}.type='uniform';
prior_all{im_all}{im}.name='sill';
prior_all{im_all}{im}.min=0.001;
prior_all{im_all}{im}.max=0.006;
prior_all{im_all}{im}.prior_master=i_master;

%% MATERN TYPE
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian - Matern';

prior_all{im_all}{im}.type='fftma';
prior_all{im_all}{im}.m0=m0;
prior_all{im_all}{im}.Cm=sprintf('.0003 Mat(%g,90,.3)',cov_range);
prior_all{im_all}{im}.fftma_options.pad_x=150; % avoid striping
prior_all{im_all}{im}.fftma_options.pad_y=150; % avoid striping
i_master=im;

% NU / Matern - variable parameter
im=im+1;
prior_all{im_all}{im}.type='uniform';
prior_all{im_all}{im}.name='nu';
prior_all{im_all}{im}.min=.1;
prior_all{im_all}{im}.max=2;
prior_all{im_all}{im}.prior_master=i_master;

% select the prior
prior=prior_all{use_prior};

%% PLOT SAMPLE FROM PRIOR_MUL
do_plot_prior_mul=0;
if do_plot_prior_mul==1;
    figure(2);clf;
    nsim=6;
    np=length(prior_all);
    for i=1:np
        prior=prior_all{i};

        for isim=1:nsim

            subplot(np,nsim,(i-1)*nsim+isim);
            [m,prior]=sippi_prior(prior);
            imagesc(prior{1}.x, prior{1}.y, m{1});
            set(gca,'FontSize',6)
            axis image;
            caxis(prior{1}.cax);
            if isim==1;
                t=title(sprintf('%s) %s',char(i+96),prior{1}.title));
                pos=get(t,'position');
                pos(1)=-5;
                set(t,'position',pos,'HorizontalAlignment','Left','FontSize',10)
            end
        end
        colorbar_shift;

    end
    print_mul(sprintf('%s_prior_reals',mfilename))
else
    sippi_plot_prior_sample(prior_all{use_prior});
end


%% SETUP THE FORWARD MODEL(S)
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S(i_use,:);
forward.receivers=D.R(i_use,:);
forward_ref=forward;

if use_forward==1;
    forward.type='ray_2d';
    forward.r=1;
    forward.name=forward.type;
elseif use_forward==2;
    forward.type='ray';
    forward.linear=1;
    forward.name='SR';
elseif use_forward==3;
    forward.type='ray';
    forward.linear=0;
    forward.name='BR';
elseif use_forward==4;
    forward.type='fat';
    forward.linear=1;
    forward.freq=0.1;
    forward.name='SF';
elseif use_forward==5;
    forward.type='fat';
    forward.linear=0;
    forward.freq=0.1;
    forward.name='BF';
elseif use_forward==6;
    forward.type='eikonal';
    forward.name=forward.type;
elseif use_forward==7;
    forward.type='fd';
    forward.name=forward.type;
end



if use_forward_ref==1;
    forward_ref.type='ray_2d';
    forward_ref.r=1;
    forward_ref.name=forward.type;
elseif use_forward_ref==2;
    forward_ref.type='ray';
    forward_ref.linear=1;
    forward_ref.name='SR';
elseif use_forward_ref==3;
    forward_ref.type='ray';
    forward_ref.linear=0;
    forward_ref.name='BR';
elseif use_forward_ref==4;
    forward_ref.type='fat';
    forward_ref.linear=1;
    forward_ref.freq=0.1;
    forward_ref.name='SF';
elseif use_forward_ref==5;
    forward_ref.type='fat';
    forward_ref.linear=0;
    forward_ref.freq=0.1;
    forward_ref.name='BF';
elseif use_forward_ref==6;
    forward_ref.type='eikonal';
    forward_ref.name=forward_ref.type;
elseif use_forward_ref==7;
    forward_ref.type='fd';
    forward_ref.name=forward.type;
end


%% SETUP DATA
D=load('AM13_data.mat');

id=1;
data{id}.d_obs=D.d_obs(i_use);
%data{id}.d_std=D.d_std(i_use);
data{id}.d_std=D.d_std(i_use)*0+0.1;
%data{id}.d_std=D.d_std(i_use)*0+0.01;
if use_correlated_noise==1
    data{id}.Ct=D.Ct; % Correlated noise model according to Cordua et al (2008; 2009)
end
options.txt=sprintf('%s_std%d',options.txt,ceil(10*data{1}.d_std(1)));

if use_reference==1;
    rng('default');
    rng(rseed);

    m_ref=sippi_prior(prior);
    [d_ref,forward_ref]=sippi_forward(m_ref,forward_ref,prior,data);
    [d,forward]=sippi_forward(m_ref,forward,prior,data);
    [logL,L,data]=sippi_likelihood(d_ref,data);

    d_obs=d_ref{1};
    d_noise = randn(size(data{1}.d_std)).*data{1}.d_std;
    data{1}.d_noisefree=d_obs;
    data{1}.d_obs=data{1}.d_noisefree+d_noise;

    options.mcmc.m_ref=m_ref;

    options.txt=sprintf('%s_ref%d',options.txt,rseed);
end

figure(1);
plot(data{1}.d_obs,'k.')
hold on
plot(d_ref{1},'r-')
plot(d{1},'b-')
hold off
xlabel('Data #')
ylabel('Travel time (mus)')
legend('d_{obs}','d_{ref}','d_{use}')

%% TEST THE SETUP

% generate a realization from the prior
[m,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m);
% Compute the forward response related to the realization of the prior model generated above
[d,forward]=sippi_forward(m,forward,prior,data);
try;
    forward.G=sparse(forward.G);
end

% plot the geomoetry
sippi_plot_traveltime_kernel(forward,prior,m);
figure(100);print_mul('AM13_SourceReceiverLoc');
figure(101);print_mul('AM13_Kernel');

% Compute the likelihood
[logL,L,data]=sippi_likelihood(d,data)

% plot the forward response and compare it to the observed data
sippi_plot_data(d,data);

%% SETUP MODELING ERRROR

%% MAKE SURE TO PROPAGATE ALSO THE MEASUREMENT NOISE TO NORMAL SCORE SPACE
% Such then when comjputeing the likelihood, is is dnoe conditoíonal to
% both measurement and modeling errors!!!
% This is importnat when setting up data{1}.dt and data{1}.Ct. such that
% the only chage to sippi_likelihood is to optionally apply the normal score
% residual data before likelihood computation!
% residual
%

%Nme=2000;
[Ct,dt,dd,d_full,d_app]=sippi_compute_modelization_forward_error(forward_ref,forward,prior,Nme,d,0);
[Ct_ns,dt_ns,dd_ns,d_full_ns,d_app_ns,o_nscore,dd_org_ns]=sippi_compute_modelization_forward_error(forward_ref,forward,prior,Nme,d,data,1);

save('-v7.3',sprintf('%s_MODELERRROR_N%d',options.txt,Nme))

%% Plot the infered Gaussian model for modeling error
close all
figure(21);set_paper;
subplot(4,2,1);
plot(dt{1})
ylim([-3 3])
ylabel('mean (d_{ref}-d)')
title('d_T (Gaussian \Theta)')
subplot(4,2,[3,5]);
imagesc(Ct{1});caxis([-1 1]);
colormap(cmap_linear)
axis image;
colorbar
title('C_T (Gaussian \Theta)')

subplot(4,2,2);
plot(dt_ns{1})
ylim([-5 5])
ylabel('mean (d_{ref}-d) - normal score')
title('d_T (Non-Gaussian \Theta)')
subplot(4,2,[4,6]);
imagesc(Ct_ns{1});caxis([-1 1]);
colormap(cmap_linear)
axis image;
colorbar
title('C_T (Non-Gaussian \Theta)')
print_mul(sprintf('%s_dtCT',options.txt))


%% Actual Modeling Error
figure(22);
j=0;
for id=ceil(linspace(1,length(d{1}),9));
    j=j+1;
    subplot(3,3,j);
    h=histogram(dd{1}(id,:),30,'FaceColor',[0 0 0]);
    BinCenter=(h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;    
    BinWidth=h.BinWidth;
    BinPDF=(1/BinWidth).*h.BinCounts./sum(h.BinCounts);
    bar(BinCenter,BinPDF,'FaceColor',[0 0 0]+.3)

    m0=dt{1}(id);
    s0=sqrt(Ct{1}(id,id));
    BinCenter_G = linspace(m0-3*s0,m0+3*s0,51);
    BinPDF_G=normpdf(BinCenter_G,m0,s0);
    hold on
    plot(BinCenter_G,BinPDF_G,'r-','LineWidth',2)
    hold off

    title(sprintf('id=%d',id))
end
sgtitle('Actual modeling error')
print_mul(sprintf('%s_MEsample',options.txt))


figure(23);
j=0;
for id=ceil(linspace(1,length(d{1}),9));
    j=j+1;
    subplot(3,3,j);
    h=histogram(dd_ns{1}(id,:),30,'FaceColor',[0 0 0]);
    BinCenter=(h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;    
    BinWidth=h.BinWidth;
    BinPDF=(1/BinWidth).*h.BinCounts./sum(h.BinCounts);
    bar(BinCenter,BinPDF,'FaceColor',[0 0 0]+.3)


    m0=dt_ns{1}(id);
    s0=sqrt(Ct_ns{1}(id,id));
    BinCenter_G = linspace(m0-3*s0,m0+3*s0,51);
    BinPDF_G=normpdf(BinCenter_G,m0,s0);
    hold on
    plot(BinCenter_G,BinPDF_G,'r-','LineWidth',2)
    hold off


    %hist(dd_ns{1}(id,:),30);
    title(sprintf('id=%d',id))
end
sgtitle('Actual modeling error in normal score space')
print_mul(sprintf('%s_MEsampleNS',options.txt))


%% SIMULATE MODELING ERRORS
nd=length(d{1});
nr = Nme;
irays=[1,50,100,250,500];
d_real = gaussian_simulation_cholesky(dt{1},Ct{1}+0.001*eye(nd),nr);
d_real_ns_ns = gaussian_simulation_cholesky(dt_ns{1},Ct_ns{1}+0.001*eye(nd),nr);
d_real_ns = inscore_mul(d_real_ns_ns,o_nscore{1});

figure(21)
subplot(1,1,1);
hold on
plot(d{1}+dd{1}(:,1:nr),'b-')
plot(d{1}+d_real(:,1:nr),'r-')
plot(d{1},'k-','LineWidth',3)
hold off
grid on
sgtitle('Real and simulated Gaussian modeling error')


figure(22);
subplot(1,1,1);
plot(d{1},'k-','LineWidth',3)
hold on
plot(d{1}+dd{1}(:,1:nr),'b-','LineWidth',.1)
plot(d{1}+d_real_ns(:,1:nr),'r-','LineWidth',.1)
plot(d{1},'k-','LineWidth',1)
hold off
grid on
legend({'d','dd real','dd sim'})
sgtitle('Real and simulated non-Gaussian modeling error')

figure(23);
j=0;
for id=ceil(linspace(1,length(d{1}),9));
    j=j+1;
    subplot(3,3,j);
    [hn,hc]=hist([dd{1}(id,:);d_real(id,:),;d_real_ns(id,:)]',30);
    p=plot(hc,hn,'-');
    set(p(1),'LineWidth',3)
    title(sprintf('id=%d',id))
    legend('True','G','NS')
end
sgtitle('Modeling error - real and simulated')


%% Compute average log-likelihood from all computes data residuals!
iCt = inv(Ct{1}+diag(data{1}.d_std.^2));
iCt_ns = inv(Ct_ns{1});
for ir=1:min([10 Nme])
    delta_d = dd{1}(:,ir);

    % Only correlated noise
    d_std = data{1}.d_std;
    logL_mul(ir,1)=sum(-.5*delta_d.^2./d_std.^2);

    % Gaussian ME
    dd0=delta_d-dt{1};
    logL_mul(ir,2)=-.5*dd0'*iCt*dd0;

    % Gaussian ME
    dd0_ns=nscore_mul(delta_d,o_nscore{1})-dt_ns{1};
    logL_mul(ir,3)=-.5*dd0_ns'*iCt_ns*dd0_ns;

    figure_focus(35);
    subplot(1,2,1);
    plot(1:ir,logL_mul(1:ir,2:3),'-');
    legend({'ME','ME_{NS}'})
    subplot(1,2,2);
    plot(1:ir,logL_mul(1:ir,:),'-');
    legend({'Normal','ME','ME_{NS}'})
    ylabel(log(L))
    drawnow;

end
print_mul(sprintf('%s_logLtrain',options.txt))

%% Test likelihood - should be really bad when using the Gaussian modeling error!!
%make sure the following forward normal score works:
% SHould we smooth the normal score a bit??

%m=sippi_prior(prior);
d_ref=sippi_forward(m_ref,forward_ref,prior);
d=sippi_forward(m_ref,forward,prior);

data_org = data;
clear data
id=1;
data{id}.d_obs=d_ref{1};
data{id}.d_std=d_ref{1}.*0+data_org{1}.d_std;

data_me = data;
data_me{1}.dt=dt{1};
data_me{1}.Ct=Ct{1};

data_ns = data;
data_ns{1}.d_std=0*data_ns{1}.d_std;
data_ns{1}.n_score=o_nscore{1};
data_ns{1}.dt=dt_ns{1};
data_ns{1}.Ct=Ct_ns{1};


tic;logL(1)=sippi_likelihood(d,data);t(1)=toc;
tic;logL(2)=sippi_likelihood(d,data_me);t(2)=toc;
tic;logL(3)=sippi_likelihood(d,data_ns);t(3)=toc;
disp(sprintf('logL=%g ',logL))

figure(13);
plot(d{1},'r-')
hold on
plot(d_ref{1},'k-')
hold off
grid on
%% TEST PROB INV (ON BOTH SYNTH AND REAL DATA) USING
% Case 1: Ignoring modeling error
% Case 2: Using Gaussian modeling error
% Case 3: Using nonGaussian modeling error

options.mcmc.nite=n_ite;
options.mcmc.i_plot=max([2500 ceil(n_ite/100)]);
options.mcmc.i_sample=options.mcmc.nite/n_reals_out;

% ANNEALING
% example of starting with a high temperature, that allow higher
% exploration. The temperature is lowered to T=1, after which the
% algorithm proceeds as a usual Metropolis sampler
if doAnneal==1;
    options.txt=[options.txt,'_','anneal'];
    i_stop_anneal=max([1000 ceil(n_ite/10)]);
    for im=1:length(prior);
        prior{im}.seq_gibbs.i_update_step_max=2*i_stop_anneal;
    end
    options.mcmc.anneal.T_begin=25; % Start temperature
    options.mcmc.anneal.T_end=1; % End temperature, at ptions.mcmc.anneal.
    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=i_stop_anneal; %  iteration number when annealing stops
end

% TEMPERING
% example of using parallel tempering (Sambridge, 2013)
if doTempering==1;
    options.txt=[options.txt,'_','temper'];
    options.mcmc.n_chains=4; % set number of chains (def=1)
    options.mcmc.T=[1 1.25 1.5 1.75]; % set number of chains (def=1)
end


O=options;
Omul{1}=sippi_metropolis(data,prior,forward,O);
O_me=options;O_me.txt=[O.txt,'_me'];
Omul{2}=sippi_metropolis(data_me,prior,forward,O_me);
O_ns=options;O_ns.txt=[O.txt,'_ns'];
Omul{3}=sippi_metropolis(data_ns,prior,forward,O_ns);

save(options.txt)

% Test likelihood of reference model versus liellhood of samples of the
% posterior.
% Hopefully ME+NS does a better job than ME alone???
% Should we use a single normal score for all travektime residualkt or (as
% now) one normal socre per data?? What goves most robist results!!
% Does the data residual using ME giev positive data residuals (it should
% not)

%%
for i=1:length(Omul);
    logL_out(i,:)=Omul{i}.C{1}.mcmc.logL;
    [reals{i},etype_mean{i},etype_var{i},reals_all{i},reals_ite{i}]=sippi_get_sample(Omul{i}.txt);
end

for i=1:length(Omul);
    for ir=1:size(reals_all{1},1)
        v = reshape(reals_all{i}(ir,:),prior{1}.dim(2),prior{1}.dim(1));
        cc=corrcoef(v(:),m_ref{1}(:));cc=cc(2);
        dv=v(:)-m_ref{1}(:);
        mae_out(i,ir)=mean(abs(dv));
        mse_out(i,ir)=sqrt(mean(dv.^2));
        cc_out(i,ir)=cc;
    end
    i1=ceil(ir/2);
    mean_cc(i)=mean(cc_out(i,i1:end))
    mean_mse(i)=mean(mse_out(i,i1:end))
    mean_mae(i)=mean(mae_out(i,i1:end))
end

figure(1);set_paper;
cax=[0.14]+[-1 1].*0.005;
for i=1:length(Omul);
    subplot(2,4,i)
    imagesc(prior{1}.x,prior{1}.y,etype_mean{i});axis image
    caxis(cax)
    if i==1;
        title('No modeling error')
    elseif i==2
        title('Gaussian modeling error')
    elseif i==3
        title('Normal Score modeling error')
    end
    subplot(2,4,i+4)
    imagesc(prior{1}.x,prior{1}.y,sqrt(etype_var{i}));axis image
    caxis([0 0.01])
    end
subplot(2,4,4)
imagesc(prior{1}.x,prior{1}.y,options.mcmc.m_ref{1});axis image
title('Reference')
caxis(cax)
print_mul(sprintf('%s_compare_mean',options.txt))

figure(33);
plot(logL_out');
legend({'Standard','ME','ME_{NS}'})
print_mul(sprintf('%s_logL',options.txt))

%%
figure(34);
ax1=subplot(3,1,1)
p=plot(cc_out');
hold on
tmp=cc_out;for i=1:size(cc_out,1);tmp(i,:)=mean_cc(i);end
p2=plot(tmp','-','LineWidth',2);
hold off
for i=1:length(p2);set(p2(i),'Color',p(i).Color);end
legend(p,{'No C_T','Gauss','Ns'});grid on
ylabel('Correlation Coefficient - to m-ref');

ax2=subplot(3,1,2)
p=plot(mae_out');
hold on
tmp=mae_out;for i=1:size(mae_out,1);tmp(i,:)=mean_mae(i);end
p2=plot(tmp','-','LineWidth',2);
for i=1:length(p2);set(p2(i),'Color',p(i).Color);end
hold off
legend(p,{'No C_T','Gauss','Ns'});grid on
ylabel('Mean Absolute Error')

ax3=subplot(3,1,3)
p=plot(mse_out');
hold on
tmp=mse_out;for i=1:size(mse_out,1);tmp(i,:)=mean_mse(i);end
p2=plot(tmp','-','LineWidth',2);
for i=1:length(p2);set(p2(i),'Color',p(i).Color);end
hold off
ylabel('Mean Squared Error')
xlabel('Realization number')
legend(p,{'No C_T','Gauss','Ns'});grid on
ylim(1)=ax2.YLim(1);
ylim(2)=ax3.YLim(2);
ax2.YLim=ylim;
ax3.YLim=ylim;
linkaxes([ax1, ax2, ax3],'x');
sgtitle(options.txt,'Interpreter','None')
print_mul(sprintf('%s_cc',options.txt))

%%

return

%% REJECTION
use_rejection=1;
if use_rejection==1;
    Nl=200000;
    Nl_me=0;
    ABC=sippi_abc_setup(prior,forward,Nl,Nl_me);

    %%
    cax=[0.12 0.17];
    cax=[0.13 0.16];
    cmap=jet;
    cax_std=[0 0.01];
    ns=100;
    j=0;
    for nd_use = [1:70]
        j=j+1;
        i_use = ceil(linspace(1,length(data{1}.d_obs),nd_use));
        %i_use = 1:nd_use;
        data_me{1}.i_use = i_use;
        data_ns{1}.i_use = i_use;
        data{1}.i_use = i_use;


        ABC.use_sippi_likelihood=0;
        data{1}.full_likelihood=0;
        data_ns{1}.full_likelihood=ABC.use_sippi_likelihood;
        data_me{1}.full_likelihood=ABC.use_sippi_likelihood;

        [logL,evidence,T_est]=sippi_abc_logl(ABC,data);
        T_est = 10;
        [m_real]=sippi_abc_post_sample(ABC,ns,T_est,logL);
        [m_mean,m_var]=sippi_abc_post_mean(ABC,T_est,logL);

        [logL_me,evidence_me,T_est_me]=sippi_abc_logl(ABC,data_me);
        T_est_me = 10;
        [m_real_me]=sippi_abc_post_sample(ABC,ns,T_est_me,logL_me);
        [m_mean_me,m_var_me]=sippi_abc_post_mean(ABC,T_est_me,logL_me);

        [logL_ns,evidence_ns,T_est_ns]=sippi_abc_logl(ABC,data_ns);
        T_est_ns = 10;
        [m_real_ns]=sippi_abc_post_sample(ABC,ns,T_est_ns,logL_ns);
        [m_mean_ns,m_var_ns]=sippi_abc_post_mean(ABC,T_est_ns,logL_ns);

        for is=1:ns
            cc=corrcoef(m_ref{1}(:),m_real{1}(:,is));cc=cc(2);
            cc_me=corrcoef(m_ref{1}(:),m_real_me{1}(:,is));cc_me=cc_me(2);
            cc_ns=corrcoef(m_ref{1}(:),m_real_ns{1}(:,is));cc_ns=cc_ns(2);

            C(is,1,j)=cc;
            C(is,2,j)=cc_me;
            C(is,3,j)=cc_ns;
        end
        meanC(j,:)= mean(C(:,:,j));




        figure(41);clf;
        subplot(2,4,1);
        imagesc(prior{1}.x,prior{1}.y,reshape(m_mean{1},prior{1}.dim(2),prior{1}.dim(1)))
        axis image;caxis(cax);colormap(gca,cmap)
        title(sprintf('No Me, T=%2.1f',T_est))
        subplot(2,4,2);
        imagesc(prior{1}.x,prior{1}.y,reshape(m_mean_me{1},prior{1}.dim(2),prior{1}.dim(1)))
        axis image;caxis(cax);colormap(gca,cmap)
        title(sprintf('Me, T=%2.1f',T_est_me))
        subplot(2,4,3);
        imagesc(prior{1}.x,prior{1}.y,reshape(m_mean_ns{1},prior{1}.dim(2),prior{1}.dim(1)))
        axis image;caxis(cax);colormap(gca,cmap)
        title(sprintf('Me-NS, T=%2.1f',T_est_ns))
        subplot(2,4,4);
        imagesc(prior{1}.x,prior{1}.y,m_ref{1})
        axis image;caxis(cax);colormap(gca,cmap)

        subplot(2,4,5);
        imagesc(prior{1}.x,prior{1}.y,reshape(sqrt(m_var{1}),prior{1}.dim(2),prior{1}.dim(1)))
        axis image;cax_std=caxis;colormap(gca,hot)
        subplot(2,4,6);
        imagesc(prior{1}.x,prior{1}.y,reshape(sqrt(m_var_me{1}),prior{1}.dim(2),prior{1}.dim(1)))
        axis image;cax_std=caxis;colormap(gca,hot)
        subplot(2,4,7);
        imagesc(prior{1}.x,prior{1}.y,reshape(sqrt(m_var_ns{1}),prior{1}.dim(2),prior{1}.dim(1)))
        axis image;caxis(cax_std);colormap(gca,hot)

        subplot(2,4,8);
        sippi_plot_traveltime_kernel(forward,prior,m_ref,0,data{1}.i_use);

        sgtitle(sprintf('%s  --  nduse= %d',options.txt,nd_use),'interpreter','none')
        drawnow;

        figure(42);
        plot(meanC(1:j,:))
        legend({'NoME','ME','NS'})
        drawnow

    end

end

%%


