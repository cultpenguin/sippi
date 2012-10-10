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
%    options.mcmc.nite [1]  : Number if iterations
%    options.mcmc.i_plot [1]: Number of iterations between updating plots
%    options.mcmc.i_sample=500: Number of iterations between saving model to disk
%
%    options.mcmc.m_init : Manually chosen starting model
%    options.mcmc.m_ref  : Reference known target model
%
%
%    options_mcmc.accept_only_improvements [0] : Optimization
% See also sippi_rejection
%
options.null='';

if ~isfield(options,'txt')
    options.txt='sippi_metropolis';
end
options.txt=sprintf('%s_metropolis_%s',datestr(now,'YYYYmmdd_HHMM'),options.txt);

try
    options.txt=sprintf('%s_%s',options.txt,forward.type);
end
try
    if forward.linear==1;
        t='lin';
    else
        t='nonlin';
    end
    options.txt=sprintf('%s_%s',options.txt,t);
end

nm=length(prior);


%% INITIALIZE ASC FILE
start_dir=pwd;
try;
    mkdir(options.txt);
    cd(options.txt);
    addpath(['..',filesep])
    
    for im=1:length(prior);
        try
            if isunix
                system(sprintf('cp ..%s%s . ',filesep,prior{im}.ti));
            else
                system(sprintf('copy ..%s%s ',filesep,prior{im}.ti));
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
rand('seed',1);
if isfield(mcmc,'m_init');
    m_init=mcmc.m_init;
    disp(sprintf('Using supplied model as starting model',mfilename))
else
    [m_init,prior] = sippi_prior(prior);
end
m_current=m_init;





%% INITIAL LIKELIHOODS
if isfield(forward,'forward_function');
    [d_init,forward,prior,data]=feval(forward.forward_function,m_init,forward,prior,data);
else
    [d_init,forward,prior,data]=sippi_forward(m_init,forward,prior,data);
end
[logL_init,L_init,data]=sippi_likelihood(d_init,data);
logL_current=logL_init;

%% COMPUTE TIME PER ITERAION

% COMPUTE THE TIME OF ONE CALL TO SIPPI_PRIOR
[m_tmp,prior_tmp] = sippi_prior(prior); % make sure prior is set
tic
[m_tmp,prior_tmp] = sippi_prior(prior_tmp);
t_prior=toc;
clear m_tmp prior_tmp;

% COMPUTE THE TIME OF ONE CALL TO SIPPI_FORWARD
tic
if isfield(forward,'forward_function');
    [d_init,forward,prior,data]=feval(forward.forward_function,m_init,forward,prior,data);
else
    [d_init,forward,prior,data]=sippi_forward(m_init,forward,prior,data);
end
t_data=toc;

% compute the number of iterations between text updates on the screen based
% on the of one forward evaluation (t_data), if not alloready set
if isfield(options,'i_update_txt');
    i_update_txt=options.i_update_txt;
else
    i_update_txt=max([1 ceil(5./(t_data+t_prior))]);
end

%% PRE ALLOCATE ARRAY FOR MCMC OUTPUT
mcmc.logL=zeros(1,mcmc.nite);
mcmc.acc=zeros(nm,mcmc.nite);
mcmc.perturb=zeros(nm,mcmc.nite);
mcmc.step=zeros(nm,mcmc.nite);

N_post_reals=floor(mcmc.nite/mcmc.i_sample);
mcmc.i_sample_logL=zeros(1,N_post_reals);

if mcmc.store_all==1
    mcmc.logL_pro = zeros(1,mcmc.nite);
    mcmc.m_pro = zeros(nm,mcmc.nite);
end


%% START THE METROPOLOS ALGORITHM
t0=now;
iacc=0;
isample=0;
for i=1:mcmc.nite;
    mcmc.i=i;
    
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
    % solve forward problem, and compute likelihood
    if isfield(forward,'forward_function');
        [d,forward,prior_propose,data]=feval(forward.forward_function,m_propose,forward,prior_propose,data);
    else
        [d,forward,prior_propose,data]=sippi_forward(m_propose,forward,prior_propose,data);
    end
    [logL_propose,L_propose,data]=sippi_likelihood(d,data);

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
    
    if Pacc>rand(1)
        % ACCEPT MODEL
        if ((i/i_update_txt)==round(i/i_update_txt))
            [t_end_txt,t_left_seconds]=time_loop_end(t0,i,mcmc.nite);
            disp(sprintf('%06d/%06d (%10s): acc %5g %5g ',mcmc.i,mcmc.nite,t_end_txt,logL_current,logL_propose))
        end
        
        prior=prior_propose; % NEEDED FOR GAUSSIAN TYPE PRIOR
        m_current=m_propose;
        logL_current=logL_propose;
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
    if ((mcmc.i/(10*mcmc.i_sample))==round( mcmc.i/(10*mcmc.i_sample) ))
        save(filename_mat)
    end
    
    %% plot current model
    if ((mcmc.i/mcmc.i_plot)==round( mcmc.i/mcmc.i_plot ))
        try
            sippi_plot_current_model(mcmc,data,d,m_current,prior);
        catch
            disp(sprintf('%s : Could not plot current model info',mfilename))
        end
        drawnow;
    end
end
save(filename_mat)
disp(sprintf('%s : DONE McMC on %s',mfilename,options.txt))

%% PLOT STATS

try;
    %sippi_plot_prior(options.txt);
catch
    disp('Could not plot figures on prior statistics')
end
try;
    %sippi_plot_posterior(options.txt);
catch
    disp('Could not plot figures on posterior statistics')
end

%% return mcmc options
options.mcmc=mcmc;

%%
cd(start_dir);

