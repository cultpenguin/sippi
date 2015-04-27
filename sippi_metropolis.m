function [options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)
% sippi_metropolis Extended Metropolis sampling in SIPPI
%
% Metropolis sampling.
%   See e.g. Hansen, T. M., Cordua, K. S., and Mosegaard, K., 2012. 
%     Inverse problems with non-trivial priors - Efficient solution through Sequential Gibbs Sampling. 
%     Computational Geosciences. doi:10.1007/s10596-011-9271-1.
%
% Call :
%    [options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)
% Input : 
%    data : sippi data structure
%    prior : sippi prior structure
%    forward : sippi forward structure
%
% options : 
%    options.txt [string] : string to be used as part of all output files
%
%    options.mcmc.nite=30000;   % [1] : Number if iterations
%    options.mcmc.i_sample=100; % : Number of iterations between saving model to disk
%    options.mcmc.i_plot=50;  % [1]: Number of iterations between updating plots
%    options.mcmc.i_save_workspace=10000;  % [1]: Number of iterations between
%                                            saving the complete workspace
%    options.mcmc.i_sample=100; % : Number of iterations between saving model to disk
%
%    options.mcmc.m_init : Manually chosen starting model
%    options.mcmc.m_ref  : Reference known target model
%
%    options_mcmc.accept_only_improvements [0] : Optimization
%
%   %% PERTUBATION STRATEGY
%   options.mcmc.pert_strategy.perturb_all=1; % Perturb all priors in each 
%                                              % iteration. def =[0]
%    %% SIMULATED ANNEALING 
%    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
%    options.mcmc.anneal.i_end=100000; %  iteration number when annealing stops
%    options.mcmc.anneal.fac_begin=20; % default, noise is scaled by fac_begin at iteration i_begin
%    options.mcmc.anneal.fac_end=1; % default, noise is scaled by fac_end at iteration i_end
%
%
% See also sippi_rejection
%
options.null='';
if ~isfield(options,'txt');options.txt='';end
if ~isempty(options.txt)
    options.txt=sprintf('%s_sippi_metropolis_%s',datestr(now,'YYYYmmdd_HHMM'),options.txt);
else
    options.txt=sprintf('%s_sippi_metropolis',datestr(now,'YYYYmmdd_HHMM'));
end
try
    options.txt=sprintf('%s_%s',options.txt,forward.type);
catch
    % No forward.type
end
try
    if forward.linear==1;
        t='lin';
    else
        t='nonlin';
    end
    options.txt=sprintf('%s_%s',options.txt,t);
catch
    % no obsolete, only applies when using the traveltime toolbox. Perhaps
    % remove?
end

nm=length(prior);


%% INITIALIZE ASC FILE
start_dir=pwd;
try;
    mkdir(options.txt);
    cd(options.txt);
    addpath(['..',filesep])
    
    for im=1:length(prior);
        if isstr(prior{im}.ti)
            try
                if isunix
                    system(sprintf('cp ..%s%s . ',filesep,prior{im}.ti));
                else
                    system(sprintf('copy ..%s%s ',filesep,prior{im}.ti));
                end
            end
        end
    end
end
for im=1:nm
    filename_asc{im}=sprintf('%s_m%d%s',options.txt,im,'.asc');
    fid=fopen(filename_asc{im},'w');
    fclose(fid);
end
filename_mat=[options.txt,'.mat'];

% SET DFAULT PLOTTING SETTINGS
options=sippi_plot_defaults(options);


%% INITIALIZE MCMC OPTIONS
options=sippi_mcmc_init(options,prior);
mcmc=options.mcmc;


%% INITIALIZE prior
prior=sippi_prior_init(prior);

%% INITIALIZE data
%data=sippi_data_init(data);

for id=1:length(data)
    if ~isfield(data{id},'i_use');
        data{id}.i_use=find(ones(size(data{id}.d_obs)));
    end
    %data{id}.N=prod(size(data{id}.d_obs));
    data{id}.N=length(data{id}.i_use);
end
%% INITIALIZE data/likelihood
% data=sippi_data_init(data)



%% STARTING  MODEL
if isfield(mcmc,'m_init');
    m_init=mcmc.m_init;
    disp(sprintf('Using supplied model as starting model',mfilename))
else
    [m_init,prior] = sippi_prior(prior);
end
m_current=m_init;


%% INITIAL LIKELIHOODS
[d_init,forward,prior,data]=sippi_forward(m_init,forward,prior,data);

%% Initialize
data_init=data;
if isfield(mcmc,'anneal');
    % SET ANNEALING IF APPLICABLE
    for id=1:length(data_init);data_init{id}.recomputeCD=1;end
    [data_init,mcmc]=sippi_anneal_adjust_noise(data_init,1,options.mcmc,prior);
end
[logL_init,L_init]=sippi_likelihood(d_init,data_init);
logL_current=logL_init;
L_current=L_init;

%% COMPUTE TIME PER ITERAION
% COMPUTE THE TIME OF ONE CALL TO SIPPI_PRIOR
[m_tmp,prior_tmp] = sippi_prior(prior); % make sure prior is set
tic
[m_tmp,prior_tmp] = sippi_prior(prior_tmp);
t_prior=toc;
clear m_tmp prior_tmp;

% COMPUTE THE TIME OF ONE CALL TO SIPPI_FORWARD
[d_init,forward,prior,data]=sippi_forward(m_init,forward,prior,data);
d_current=d_init;
t_data=toc;

% compute the number of iterations between text updates on the screen based
% on the of one forward evaluation (t_data), if not alloready set
if isfield(options,'i_update_txt');
    i_update_txt=options.i_update_txt;
else
    i_update_txt=max([1 ceil(10./(t_data+t_prior))]);
end

%% PRE ALLOCATE ARRAY FOR MCMC OUTPUT
mcmc.logL=zeros(1,mcmc.nite);
if length(data)>1
    mcmc.logL_all=zeros(length(data),mcmc.nite);
end
mcmc.acc=zeros(nm,mcmc.nite);
mcmc.perturb=zeros(nm,mcmc.nite);
mcmc.step=zeros(nm,mcmc.nite);
mcmc.time=zeros(1,mcmc.nite);

N_post_reals=floor(mcmc.nite/mcmc.i_sample);
mcmc.i_sample_logL=zeros(1,N_post_reals);

if mcmc.store_all==1
    mcmc.logL_pro = zeros(1,mcmc.nite);
    mcmc.m_pro = zeros(nm,mcmc.nite);
end

%% START THE METROPOLOS ALGORITHM
disp(sprintf('%s : starting extended Metropolis sampler in %s',mfilename,options.txt))
t0=now;
iacc=0;
isample=0;
for i=1:mcmc.nite;
    mcmc.i=i;
    mcmc.time(i)=now;
    % set seed
    for im=1:length(prior)
        prior{im}.seed=i;
    end
    
    %% Next section only for keeping track of FFT-MA options / Seg Gibbs
    for im=1:length(prior);
        if isfield(prior{im},'fftma_options');
            fftma_options{im}=prior{im}.fftma_options;
        end
    end

    %% SELECT IF STEP LENGTH HAS TO BE UPDATED
    for im=1:length(prior)
        if (((mcmc.i./prior{im}.seq_gibbs.i_update_step)==round((mcmc.i./prior{im}.seq_gibbs.i_update_step)))&(mcmc.i<prior{im}.seq_gibbs.i_update_step_max))
            % UPDATE STEP LENGTH
            prior=sippi_prior_set_steplength(prior,mcmc,im);
        end
        mcmc.step(im,i)=prior{im}.seq_gibbs.step(1);
    end
    
    %% Sample prior
    % SELECT WHICH MODEL PARAMETERS TO PERTURB

    for im=1:length(prior);
        prior{im}.perturb=0;
    end
    
    % PERTURBATION STRATEGY
    if mcmc.pert_strategy.perturb_all==1,
        im_perturb=1:1:length(prior);
    else
        % perturb one parameter according to frequency distribuiton
        i_pert=mcmc.pert_strategy.i_pert;
        pert_freq=cumsum(mcmc.pert_strategy.i_pert_freq);
        pert_freq=pert_freq./max(pert_freq);
        im_perturb=i_pert(min(find(rand(1)<pert_freq)));
    end
   
    for k=1:length(im_perturb);
        prior{im_perturb(k)}.perturb=1;
        mcmc.perturb(im_perturb(k),mcmc.i)=1;
    end
    
    
    % Possibly pertrub more than one model parameter
    %prior{im_perturb}.perturb=1;
    %mcmc.perturb(im_perturb,mcmc.i)=1;
 
    % SAMPLE PRIOR
    [m_propose,prior_propose] = sippi_prior(prior,m_current);
    
    %% FORWARD PROBLEM
    [d,forward,prior_propose,data]=sippi_forward(m_propose,forward,prior_propose,data);
    do_anneal=0;
    if isfield(mcmc,'anneal');
        if (i>=mcmc.anneal.i_begin)&(i<=mcmc.anneal.i_end)
            do_anneal=1;
        end
    end 
    if do_anneal==1,
        % DO ANNEAL
        [data_test]=sippi_anneal_adjust_noise(data,i,options.mcmc,prior);
        for id=1:length(data_test);data_test{id}.recomputeCD=1;end
        [logL_propose,L_propose]=sippi_likelihood(d,data_test);
        [logL_current,L_current]=sippi_likelihood(d_current,data_test);
    else
        [logL_propose,L_propose,data]=sippi_likelihood(d,data);
    end
    
    if mcmc.store_all==1
        mcmc.logL_pro(i)=logL_propose;
        for im=1:length(prior);
            mcmc.m_pro(im,i) = m_propose{im}(1);
        end
    end
    
    
    % Accept probability
    Pacc = exp(logL_propose-logL_current);
    if (mcmc.accept_only_improvements==1)
        % Optimization only?
        Pacc=(Pacc>=1); %% Pacc = 1 if Pacc>=1, else Pacc=0;
    end
    
    % Optionally accept all proposed models
    if (mcmc.accept_all==1), Pacc=1; end
    
    forward.last_proposed_model_accept=0; % 
    if Pacc>rand(1)
        % ACCEPT MODEL
        forward.last_proposed_model_accept=1;
        
        if ((i/i_update_txt)==round(i/i_update_txt))
            [t_end_txt,t_left_seconds]=time_loop_end(t0,i,mcmc.nite);
            disp(sprintf('%06d/%06d (%10s): acc %5g %5g ',mcmc.i,mcmc.nite,t_end_txt,logL_current,logL_propose))
        end
        
        prior=prior_propose; % NEEDED FOR GAUSSIAN TYPE PRIOR
        m_current=m_propose;
        d_current=d;
        logL_current=logL_propose;
        L_current=L_propose;
        iacc=iacc+1;
        mcmc.acc_logL(iacc)=logL_current;
        mcmc.acc(im_perturb,mcmc.i)=1;
        
    else
        % REJECT MODEL
        %% Next section only for keeping track of FFT-MA options / Seg Gibbs
        %  reser z_rand values
        for im=1:length(prior);            
            if isfield(prior{im},'fftma_options');prior{im}.fftma_options=fftma_options{im};end
        end
        
        %%
        if ((i/i_update_txt)==round(i/i_update_txt))
            [t_end_txt,t_left_seconds]=time_loop_end(t0,i,mcmc.nite);
            disp(sprintf('%06d/%06d (%10s):     %5g %5g %s',mcmc.i,mcmc.nite,t_end_txt,logL_current,logL_propose))
        end
    end
    mcmc.logL(i)=logL_current;
    if length(data)>1
        % store logL for all data types seperately
        mcmc.logL_all(:,i)=L_current(:);
    end
    
    % SAVE CURRENT MODEL
    if ((mcmc.i/mcmc.i_sample)==round( mcmc.i/mcmc.i_sample ))
        isample=isample+1;
        mcmc.i_sample_logL(isample)=logL_current;
        for im=1:nm
            fid=fopen(filename_asc{im},'a+');
            fprintf(fid,' %10.7g ',[m_current{im}(:)]);
            fprintf(fid,'\n');
            fclose(fid);
        end
    end
    % SAVE WORKSPACE
    if ((mcmc.i/(mcmc.i_save_workspace))==round( mcmc.i/(mcmc.i_save_workspace) ))
        try
            save(filename_mat)
        catch
            disp(sprintf('%s : failed to save data to %s',mfilename,filename_mat))
        end
    end
    
    %% plot current model
    if ((mcmc.i/mcmc.i_plot)==round( mcmc.i/mcmc.i_plot ))
        try
            sippi_plot_current_model(mcmc,data,d_current,m_current,prior);
        catch
            disp(sprintf('%s : Could not plot current model info',mfilename))
        end
        try;figure_focus(3);print_muk('current_status',options.plot.hardcopy_types);end
        drawnow;
    end
end


%% return mcmc options
mcmc.m_current=m_current;
options.mcmc=mcmc;

save(filename_mat)
disp(sprintf('%s : DONE McMC on %s',mfilename,options.txt))

%%
cd(start_dir);

