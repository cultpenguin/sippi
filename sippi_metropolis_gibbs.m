function [options,data,prior,forward,m_current]=sippi_metropolis_gibbs(data,prior,forward,options)
% sippi_metropolis Extended Metropolis sampling in SIPPI
%
% Metropolis sampling.
%   See e.g.
%     Hansen, T. M., Cordua, K. S., and Mosegaard, K., 2012.
%     Inverse problems with non-trivial priors - Efficient solution through Sequential Gibbs Sampling.
%     Computational Geosciences. doi:10.1007/s10596-011-9271-1.
%
%     Sambridge, M., 2013 - A Parallel Tempering algorithm for
%     probabilistic sampling and multi-modal optimization.
%
% Call :
%    [options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)
% Input :
%    data : sippi data structure
%    prior : sippi prior structure
%    forward : sippi forward structure
%
% options :
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
%    options.mcmc.accept_only_improvements [0] : Optimization
%    options.mcmc.accept_all [0]: accepts all proposed models (ignores lilkelihood)
%
%    options.txt [string] : string to be used as part of all output file names
%
%    %% PERTURBATION STRATEGY
%    % Perturb all model parameter all the time
%    options.mcmc.pert_strategy.perturb_all=1; % Perturb all priors in each
%                                              % iteration. def =[0]
%    options.mcmc.pert_strategy.perturb_all=2; % Perturb a random selection of
%                                              % all priors in each iteration. def =[0]
%
%    % Perturb one a prior type at a time, according to some frequency
%    options.mcmc.pert_strategy.i_pert = [1,3]; % only perturb prior 1 and 3
%    options.mcmc.pert_strategy.i_pert_freq = [2 8]; % perturb prior 3 80% of
%                                               % the time and prior 1 20%
%                                               % of the time
%    % the default pertubation strategt is to select one prior model to
%    % perturb at random for each iteration
%
%
%    %% TEMPERING
%    options.mcmc.n_chains=3; % set number of chains (def=1)
%    options.mcmc.T=[1 1.1 1.2];      % set temperature of chains [1:n_chains]
%    options.mcmc.chain_frequency_jump=0.1; % probability allowing a jump
%                                           %  between two chains
%    %% ANNEALING (TEMPERATURE AS A FUNCTION OF ITERATION NUMBER)
%    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
%    options.mcmc.anneal.i_end=100000; %  iteration number when annealing stops
%    options.mcmc.anneal.T_begin=5; % Start temperature for annealing
%    options.mcmc.anneal.T_end=1; % End temperature for annealing
%
%    %% GIBBS SAMPLING of 1D prior TYPES AT SOME ITERATIONS
%    options.mcmc.gibbs.i_gibbs = 10; % Use Gibbs sampling for every i_gibbs iterations
%                                       [-1]: no Gibbs sampling
%    options.mcmc.gibbs.usedim=2; % determines whether Gibbs sampling is 
%                                   performed on 1D[1] or 2D[2] conditional 
%    options.mcmc.gibbs.i_pert = [1 2 3]% select the prior ids to use for Gibbs sampling (must be 1D)
%                                       % if not set all 1D priors are
%                                       % considered
%    options.mcmc.gibbs.Nm=31; % number random realizaitions from the prior
%                              % used to estimate the 1D/2D marginal
%
%    %% VERBOSITY
%    The amount of text info displayed at the prompt, can be controlled by
%    setenv('SIPPI_VERBOSE_LEVEL','2') % all: information on chain swapping
%    setenv('SIPPI_VERBOSE_LEVEL','1') % information about seq-gibbs step update
%    setenv('SIPPI_VERBOSE_LEVEL','0'); % [def] frequent update
%    setenv('SIPPI_VERBOSE_LEVEL','-1'); % rare update om finish time
%    setenv('SIPPI_VERBOSE_LEVEL','-2'); % indication of stop and start
%    setenv('SIPPI_VERBOSE_LEVEL','-3'); % none
%
%    %% MULTIPLE RUNS (IN PARALLEL)
%    % In case the matlab parallel toolbox is installed, then a selected
%    % number of indenpendent]] Metropolis chains can be run in parallel
%    % using (see also sippi_metropolis_mulrun.m):
%    options.nruns = 4; % to run 4 independent runs
%
%    %% STARTING FROM A SAVED STATE
%    % If for some reason sampling is stopped before end of simulation
%    % simulation can be started again from the last saved state of the 
%    % workspacem as set in "options.mcmc.i_save_workspace=10000;"
%    % Go to the folder with the mat-file and run "sippi_metropolis"
%    % with no input arguments, to restart sampling
%
% See also sippi_metropolis_mulrun, sippi_rejection
%
%


%% CHECK FOR NO INPUT ARGUMENTS - THEN START FROM MATFILE
start_from_mat_file=0;
if nargin==0;
    sippi_verbose(sprintf('%s: no input arguments, trying to continue sampling',mfilename))
    [p,f]=fileparts(pwd);
    mat_file=[f,filesep,'mat'];
    if exist(mat_file,'file');
        disp('OK')
    else
        d=dir('*.mat');
        if length(d)==0;
            sippi_verbose(sprintf('%s: NO matfile in folder - quitting',mfilename))
            options=[];
            data=[];
            prior=[];
            forward=[];
            m_current=[];
            return
        else
            mat_file=d(1).name;
        end
    end
    
    sippi_verbose(sprintf('%s: trying to load ''%s'' and continue sampling.',mfilename,mat_file))
    try
        load(mat_file,'prior','data','forward','options');
        mcmc=options.mcmc;
        start_from_mat_file=1;
    catch
        sippi_verbose(sprintf('%s: FAILED to data form load ''%s'' and continue sampling. QUITTING',mfilename,mat_file))
        options=[];
        data=[];
        prior=[];
        forward=[];
        m_current=[];
        return
    end
end

options.null='';
%% MULTIPLE RUNS IN PARALLEL
if isfield(options,'nruns');
    if options.nruns>1
        [o_all,data,prior,forward]=sippi_metropolis_mulrun(data,prior,forward,options);
        m_current=[];
        options=o_all;
        return
    end
end

if start_from_mat_file==0;
    % ONLY DO THIS IF STARTING FROM SCRATCH
    if ~isfield(options,'txt');options.txt='';end
    if ~isempty(options.txt)
        options.txt=sprintf('sm_%s_%s',datestr(now,'YYYYmmdd_HHMM'),options.txt);
    else
        options.txt=sprintf('sm_%s',datestr(now,'YYYYmmdd_HHMM'));
    end
    
    start_dir=pwd;
    %% MAKE OUTPUT DIR
    try;
        mkdir(options.txt);
        cd(options.txt);
        addpath(['..',filesep])
        
        % copy training image file if it exists
        for im=1:length(prior);
            if ischar(prior{im}.ti)
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
    nm=length(prior); % numer prior types
    
    %% INITIALIZE MCMC OPTIONS
    options=sippi_mcmc_init(options,prior);
    mcmc=options.mcmc;
    
    %% INITIALIZE prior
    prior=sippi_prior_init(prior);
    
    %% INITIALIZE data
    for id=1:length(data)
        if ~isfield(data{id},'i_use');
            data{id}.i_use=find(ones(size(data{id}.d_obs)));
        end
        %data{id}.N=prod(size(data{id}.d_obs));
        data{id}.N=length(data{id}.i_use);
    end
    
    
    %% INITIALIZE CHAINS
    NC=options.mcmc.n_chains;
    for ic=1:NC
        if isfield(options,'prior_mul');
            sippi_verbose(sprintf('%s: Using mulitple priors from options.prior_chain{%d}',mfilename,ic));
            C{ic}.prior=options.prior_mul{ic};
        else
            C{ic}.prior=prior;
        end
        C{ic}.data=data;
        C{ic}.forward=forward;
        C{ic}.i_chain=ic;
    end
    
    for ic=1:NC
        C{ic}.T=mcmc.T(ic);
    end
    
    %% STARTING  MODEL
    for ic=1:NC
        if isfield(mcmc,'m_init');
            C{ic}.m_current=mcmc.m_init;
            sippi_verbose(sprintf('%s: Using supplied model as starting model',mfilename));
        else
            [C{ic}.m_current,C{ic}.prior] = sippi_prior(C{ic}.prior);
        end
    end
    
    %% REFERENCE LIKELIHOOD?
    if isfield(options.mcmc,'m_ref');
        try
            options.mcmc.d_ref=sippi_forward(options.mcmc.m_ref,forward,prior,data);
            options.mcmc.logL_ref=sippi_likelihood(options.mcmc.d_ref,data);
        end
    end
    
    %% INITIAL LIKELIHOODS
    for ic=1:NC
        [C{ic}.d_current,C{ic}.forward,C{ic}.prior,C{ic}.data]=sippi_forward(C{ic}.m_current,C{ic}.forward,C{ic}.prior,C{ic}.data);
        C{ic}.data_current=C{ic}.data;
        [C{ic}.logL_current,C{ic}.L_current]=sippi_likelihood(C{ic}.d_current,C{ic}.data_current);
        C{ic}.prior_current=C{ic}.prior;
    end
    
    %% COMPUTE TIME PER ITERATION
    % COMPUTE THE TIME OF ONE CALL TO SIPPI_PRIOR
    for ic=1:NC
        [m_tmp,prior_tmp] = sippi_prior(C{ic}.prior); % make sure prior is set
        tic
        [m_tmp,prior_tmp] = sippi_prior(prior_tmp);
        t_prior(ic)=toc;
        
        % COMPUTE THE TIME OF ONE CALL TO SIPPI_FORWARD
        tic
        [d_init,forward,prior,data]=sippi_forward(m_tmp,C{ic}.forward,prior_tmp,C{ic}.data);
        d_current=d_init;
        t_data(ic)=toc;
        clear m_tmp prior_tmp d_init;
    end
    t_per_ite=sum(t_data)+sum(t_prior);
    
    % compute the number of iterations between text updates on the screen based
    % on the of one forward evaluation (t_data), if not alloready set
    if isfield(options,'i_update_txt');
        i_update_txt=options.i_update_txt;
    else
        vlevel=sippi_verbose;
        if vlevel<0,
            % less txt updates when vlevel<0
            t_up=100;
        else
            % more txt updates when vlevel>0
            t_up=20;
        end
        i_update_txt=max([1 ceil(t_up./(t_per_ite))]);
    end
    
    %% PRE ALLOCATE ARRAY FOR MCMC OUTPUT
    if NC>1
        mcmc.i_swap=ones(2,mcmc.nite).*NaN;
        mcmc.n_swap=0;
        mcmc.i_chain=zeros(NC,mcmc.nite);
    end
    
    for ic=1:NC
        if length(data)>1
            C{ic}.mcmc.logL_all=zeros(length(C{ic}.data),mcmc.nite);
        end
        C{ic}.mcmc.logL=zeros(1,mcmc.nite);
        C{ic}.mcmc.acc=zeros(nm,mcmc.nite);
        C{ic}.mcmc.perturb=zeros(nm,mcmc.nite);
        C{ic}.mcmc.step=zeros(nm,mcmc.nite);
        C{ic}.mcmc.time=zeros(1,mcmc.nite);
        
        if mcmc.store_all==1
            %    C{ic}.logL_all = zeros(1,mcmc.nite);
            C{ic}.m_cur_all = zeros(nm,mcmc.nite);
        end
        
        C{ic}.iacc=0;
        C{ic}.isample=0;
        
    end
    
    %% INITIALIZE ASC FILE
    for ic=1:NC
        for im=1:nm
            C{ic}.filename_asc{im}=sprintf('%s_m%d_C%d%s',options.txt,im,ic,'.asc');
            C{ic}.fid=fopen(C{ic}.filename_asc{im},'w');
            fclose(C{ic}.fid);
        end
    end
    filename_mat=[options.txt,'.mat'];
    
    %% SET DFAULT PLOTTING SETTINGS
    options=sippi_plot_defaults(options);
    
    %% START THE METROPOLIS ALGORITHM
    sippi_verbose(sprintf('%s: Starting extended Metropolis sampler in %s',mfilename,options.txt),-2);
    if NC>1
        T_str=sprintf('%3.1f ',mcmc.T(1:NC));
        sippi_verbose(sprintf('%s: Using %d chains at temperatures: %s',mfilename,NC,T_str),-2);
    end
    mcmc.t_start=now;
    
end

i=0;
%for i=1:mcmc.nite;
while i<=mcmc.nite;
    i=i+1;
    
    mcmc.i=i;
    mcmc.time(i)=now;
    
    %% OPTIONALLY LOAD FROM MAT_FILE
    if (i==1)&&(start_from_mat_file==1);
        % LOAD STATE FROM MATLAB
        load(mat_file);
    end
    
    %% SEELECT SAMPLING METHOD
    if (mcmc.gibbs.i_gibbs>0)
        useMetropolis = min([mod(i,mcmc.gibbs.i_gibbs),1]);
    else
        % always use Metropolos
        useMetropolis = 1;
    end    
    useGibbs=1-useMetropolis;
    
    
    
    % Use extenden Metropolis sampler
    if useMetropolis        
        [C,mcmc]=sippi_metropolis_iteration(C,mcmc,i);
    end
    if useGibbs
        if mcmc.gibbs.usedim==1
            % 1D GIBBS SAMPLING
            [C,mcmc]=sippi_metropolis_gibbs_random_iteration(C,mcmc,i);
        else
            % 2D GIBBS SAMPLING
            %[C,mcmc]=sippi_metropolis_gibbs_random_iteration_2d(C,mcmc,i,1); % plot every time            
            [C,mcmc]=sippi_metropolis_gibbs_random_iteration_2d(C,mcmc,i);            
        end
    end
    
    % SAVE CURRENT MODEL
    if ((mcmc.i/mcmc.i_sample)==round( mcmc.i/mcmc.i_sample ))
        for ic=1:length(C);
        C{ic}.isample=C{ic}.isample+1;
        for im=1:length(C{ic}.m_current)
            fid=fopen(C{ic}.filename_asc{im},'a+');
            fprintf(C{ic}.fid,' %10.7g ',[C{ic}.m_current{im}(:)]);
            fprintf(C{ic}.fid,'\n');
            fclose(C{ic}.fid);
        end
        end
    end
    
    %% SAVE WORKSPACE
    if ((mcmc.i/(mcmc.i_save_workspace))==round( mcmc.i/(mcmc.i_save_workspace) ))
        try
            sippi_verbose(sprintf('%s: i=%d, saving workspace to ''%s''',mfilename,i,filename_mat))
            save(filename_mat,'-v7.3')
        catch
            sippi_verbose(sprintf('%s: failed to save data to %s',mfilename,filename_mat))
        end
    end
    
    %% DISPLAY PROGRESS AND TIME TO FINISH
    if ((i/i_update_txt)==round(i/i_update_txt))
        [t_end_txt,t_left_seconds]=time_loop_end(mcmc.t_start,i,mcmc.nite);        
        vlevel=sippi_verbose;
        if vlevel>0, NC_end=NC; else NC_end=1; end
        for ic=1:NC_end
            %try
            T = C{ic}.T*C{ic}.T_fac;
            %catch
            %    T = C{ic}.T*mcmc.T_fac;
            %end
            txt=sprintf('%06d/%06d (%10s): C%02d logL_c=%5.2f(%5.2f), T=%5.2f',mcmc.i,mcmc.nite,t_end_txt,ic,C{ic}.logL_current,C{ic}.logL_propose,T);
            sippi_verbose(sprintf('%s: %s',mfilename,txt),-1);
            % MORE information at higher verbose level
            vlevel_pacc=1;
            if vlevel>=vlevel_pacc
                for im=1:length(prior);
                    i_perturb=find(C{ic}.mcmc.perturb(im,:));
                    N=min([i 100]);
                    [P_cur, N_acc, N] = sippi_compute_acceptance_rate(C{ic}.mcmc.acc(im,i_perturb),C{ic}.prior{im}.seq_gibbs.n_update_history);
                    %[P_cur, N_acc, N] = sippi_compute_acceptance_rate(C{ic}.mcmc.acc(im,i_perturb),N);
                    txt2=sprintf(' -- im=%d, P_Acc=%3.1f%% (%d/%d)',im, 100*P_cur,N_acc,N);
                    sippi_verbose(sprintf('%s: %s',mfilename,txt2),vlevel_pacc);
                end
            end
        end
    end
    
    %% PLOT CURRENT MODEL AND STATUS
    if ((mcmc.i/mcmc.i_plot)==round( mcmc.i/mcmc.i_plot ))
        try
            C{1}.mcmc.i=mcmc.i;
            sippi_plot_current_model(C{1}.mcmc,C{1}.data,C{1}.d_current,C{1}.m_current,C{1}.prior_current,options);
        catch
            sippi_verbose(sprintf('%s: Could not plot current model info',mfilename),0);
        end
        %%
        if NC>1
            figure_focus(35);clf;
            ylim=[min(C{1}.mcmc.logL(ceil(i*.1):i)),max(C{1}.mcmc.logL(ceil(i*.1):i))];
            for ic=1:NC;
                L{ic}=sprintf('T=%3.1f',C{ic}.T);
                
                lmin=min(C{ic}.mcmc.logL(ceil(i*.1):i));
                lmax=max(C{ic}.mcmc.logL(ceil(i*.1):i));
                if lmin<ylim(1),ylim(1)=lmin;end
                if lmax>ylim(2),ylim(2)=lmax;end
                
                plot(1:i,C{ic}.mcmc.logL(1:i),'-');
                hold on
            end
            hold off
            
            try
                set(gca,'ylim',ylim)
            catch
                %% Happen in vey rary cases..
            end
            legend(L,'location','northeastoutside')
            xlabel('Iteration number')
            ylabel('log(L)')
            grid on
            
        end
        drawnow
    end
end

if NC>1
    mcmc.i_swap = mcmc.i_swap(:,1:mcmc.n_swap);
end

mcmc.t_end=now;
mcmc.time_elapsed_in_seconds=3600*24*(mcmc.t_end-mcmc.t_start);

m_current=C{1}.m_current;
mcmc.m_current=m_current;

options.C=C; % PERHAPS TOO MEMORY INTENSIVE
options.mcmc=mcmc; % PERHAPS TOO MEMORY INTENSIVE

if isoctave
    save(filename_mat)
else
    save(filename_mat,'-v7.3')
end

sippi_verbose(sprintf('%s: DONE McMC in %5.2f hours (%g minutes), %s',mfilename,mcmc.time_elapsed_in_seconds/3600,mcmc.time_elapsed_in_seconds/60,options.txt),-2);

%%
cd(start_dir);
