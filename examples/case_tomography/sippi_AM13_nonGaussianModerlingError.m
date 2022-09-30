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


%% Make some choices
if ~exist('use_prior','var')
    % use_prior=1; % Gaussian
    use_prior=2; % Gaussian with bimodal target distribution
    % use_prior=3; % Gaussian with uniform target distribution
    % use_prior=4; % Plurigaussian
    % use_prior=5; % Gaussian with variable covariance parameters
    % use_prior=6; % Matern type covariance with varying nu parameter
end
if ~exist('use_forward','var')
    use_forward=1; % ray_2d - linear straight ray
    % use_forward=2; % linear ray (using eikonal)
    % use_forward=3; % bended ray
    % use_forward=4; % linear fat
    % use_forward=5; % bended fat
    % use_forward=6; % eikonal
    % use_forward=7; % waveform FD + first arriavle
end

if ~exist('use_forward_ref','var')
    % use_forward_ref=1; % ray_2d - linear straight ray
    % use_forward_ref=2; % linear ray (using eikonal)
    % use_forward_ref=3; % bended ray
    % use_forward_ref=4; % linear fat
    % use_forward_ref=5; % bended fat
    use_forward_ref=6; % eikonal
    % use_forward_ref=7; % waveform FD + first arriavle
end


if ~exist('use_metropolis','var')
    use_metropolis=1;
end
if ~exist('use_rejection','var')
    use_rejection=0;
end

if ~exist('use_reference','var')
    use_reference=0;
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

D=load('AM13_data.mat');
options.txt='AM13_gaussian';

%% SETUP DATAwji
D=load('AM13_data.mat');

id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std.*0+0.1;
if use_correlated_noise==1
    data{id}.Ct=D.Ct; % Correlated noise model according to Cordua et al (2008; 2009)
end

figure(1);
plot(data{1}.d_obs)
xlabel('Data #')
ylabel('Travel time (mus)')

%% SETUP DIFFERENT PRIOR STRUCTURES
% define some standard values
m0=0.145;
Va=sprintf('.0003 Sph(%g,90,.3)',cov_range);

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
forward.sources=D.S;
forward.receivers=D.R;
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

forward_ref.type='eikonal';
forward_ref.name=forward.type;
if use_forward_ref==7;
    forward_ref.type='fd';
    forward_ref.name=forward.type;
end

%% TEST THE SETUP
prior=prior_all{use_prior};
% generate a realization from the prior
[m,prior]=sippi_prior(prior);
sippi_plot_model(prior,m);
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

N=2000;
[Ct,dt,dd,d_full,d_app]=sippi_compute_modelization_forward_error(forward_ref,forward,prior,N,d,0);
[Ct_ns,dt_ns,dd_ns,d_full_ns,d_app_ns,o_nscore,dd_org_ns]=sippi_compute_modelization_forward_error(forward_ref,forward,prior,N,d,data,1);


save -v7.3 MODELERRROR
%% Plot the infered Gaussian model for modeling error
close all
figure(1);
subplot(4,2,1);
plot(dt{1})
ylim([-3 3])
ylabel('mean (d_{ref}-d)')
subplot(4,2,[3,5]);
imagesc(Ct{1});caxis([-1 1]);
colormap(cmap_linear)
axis image;
colorbar

subplot(4,2,2);
plot(dt_ns{1})
ylim([-3 3])
ylabel('mean (d_{ref}-d) - normal score')
subplot(4,2,[4,6]);
imagesc(Ct_ns{1});caxis([-1 1]);
colormap(cmap_linear)
axis image;
colorbar

figure(2);
j=0;
for id=ceil(linspace(1,length(d{1}),9));
    j=j+1;
    subplot(3,3,j);
    hist(dd{1}(id,:),30);
    title(sprintf('id=%d',id))
end
sgtitle('Modeling error')


figure(3);
j=0;
for id=ceil(linspace(1,length(d{1}),9));
    j=j+1;
    subplot(3,3,j);
    hist(dd_ns{1}(id,:),30);
    title(sprintf('id=%d',id))
end
sgtitle('Modeling error normal score')

%% SIMULATE MODELING ERRORS
nd=length(d{1});
nr = N;
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


%% Test likelihood - should be really bad when using the Gaussian modeling error!!
%make sure the following forward normal score works:
% SHould we smooth the normal score a bit??


m=sippi_prior(prior);
d_ref=sippi_forward(m,forward_ref,prior);
d=sippi_forward(m,forward,prior);


clear data
id=1;
data{id}.d_obs=d_ref{1};
data{id}.d_std=d_ref{1}.*0+0.0001;

data_ns = data;
data_ns{1}.d_std=0*data_ns{1}.d_std;
data_ns{1}.n_score=o_nscore{1};
data_ns{1}.dt=dt_ns{1};
data_ns{1}.Ct=Ct_ns{1};

data_me = data;
data_me{1}.dt=dt{1};
data_me{1}.Ct=Ct{1};

sippi_likelihood(d,data_me)
sippi_likelihood(d,data_ns)

%% TEST PROB INV (ON BOTH SYNTH AND REAL DATA) USING
% Case 1: Ignoring modeling error
% Case 2: Using Gaussian modeling error
% Case 3: Using nonGaussian modeling error


    
options.mcmc.nite=n_ite;
options.mcmc.i_plot=max([500 ceil(n_ite/100)]);
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



O=sippi_metropolis(data,prior,forward,options);
O_me=sippi_metropolis(data_me,prior,forward,options);
%O_ns=sippi_metropolis(data_ns,prior,forward,options);

%options.mcmc.time_elapsed_in_seconds



return
%% SET NAME OF SIMULATION
options.txt=sprintf('AM13_I%d_CO%d_r%g',use_prior,use_correlated_noise,cov_range);
options.txt=sprintf('%s_DX%d',options.txt,1000*dx);
options.txt=[options.txt,'_',forward.name];

x=prior{1}.x;
y=prior{1}.y;
try
    G=forward.G;
try
    Cd=data{1}.CD;
    save(options.txt,'x','y','m','d','logL','G','Cd')
catch
    d_std=data{1}.d_std;
    save(options.txt,'x','y','m','d','logL','G','d_std')
end
end

%% MAKE REFERENCE MODEL
if use_reference==1;
    rng(1);
    m_ref=sippi_prior(prior);
    [d,forward]=sippi_forward(m_ref,forward,prior,data);
    [logL,L,data]=sippi_likelihood(d,data);
    
    d_obs=d{1};
    d_noise=gaussian_simulation_cholesky(0,data{1}.CD,1);
    data{1}.d_obs=d_obs+d_noise;
    
    options.mcmc.m_ref=m_ref;    
    
    options.txt=[options.txt,'_ref'];    
end

%% SETUP METROPOLIS
if use_metropolis==1
    %rng('default');
    %rng(rseed);
    
    options.mcmc.nite=n_ite;
    options.mcmc.i_plot=max([500 ceil(n_ite/100)]);
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
    
    options=sippi_metropolis(data,prior,forward,options);
    options.mcmc.time_elapsed_in_seconds
    
    % PLOT SAMPLE FROM PRIOR
    sippi_plot_prior_sample(options.txt);
    % PLOT SAMPLE AND STATS FROM POSTERIOR
    if plot_posterior_sample==1;  sippi_plot_posterior_sample(options.txt);end
    if plot_posterior_data==1;    sippi_plot_posterior_data(options.txt);end
    if plot_posterior_loglikelihood==1;    sippi_plot_posterior_loglikelihood(options.txt);end
    if plot_posterior_2d_marg==1;    sippi_plot_posterior_2d_marg(options.txt);end

    
    %% PLOT PRIOR AND POSTERIO MOVIE
    % sippi_plot_movie(options.txt)
    
    
end
return

%% REJECTION
if use_rejection==1;
    rng('default');
    rng(1);
    options.mcmc.nite=n_ite;
    options.mcmc.i_plot=max([500 ceil(n_ite/100)]);
    
    options=sippi_metropolis(data,prior,forward,options);
    
    
    % PLOT SAMPLE FROM PRIOR
    sippi_plot_prior_sample(options.txt);
    % PLOT SAMPLE AND STATS FROM POSTERIOR
    if ~exist('plot_posterior_sample','var');  sippi_plot_posterior_sample(options.txt);end
    if ~exist('plot_posterior_data','var');    sippi_plot_posterior_data(options.txt);end
    if ~exist('plot_posterior_loglikelihood','var');    sippi_plot_posterior_loglikelihood(options.txt);end
    if ~exist('plot_posterior_2d_marg','var');    sippi_plot_posterior_2d_marg(options.txt);end
    
    
    
end
