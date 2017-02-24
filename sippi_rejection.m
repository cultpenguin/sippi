function options=sippi_rejection(data,prior,forward,options)
% sippi_rejection Rejection sampling
%
% Call :
%     options=sippi_rejection(data,prior,forward,options)
%
% input arguments
%
%   options.mcmc.i_plot
%   options.mcmc.nite     % maximum number of iterations
%   options.mcmc.logLmax [def=1]; % Maximum possible log-likelihood value
%                                 % (used for normaliztion)
%
%   options.mcmc.adaptive_rejection=1, adaptive setting of maximum likelihood
%                  (def=[0])
%                  At each iteration logLmax will be set if log(L(m_cur)=>options.mcmc.logLmax
%
%
%   options.mcmc.max_run_time_hours = 1; % maximum runtime in hours
%                                        % (overrides options.mcmc.nite if needed)
%
%   options.mcmc.T = 1; % Tempering temperature. T=1, implies no tempering
%
%   %% STARTING FROM A SAVED STATE
%   % If for some reason sampling is stopped before end of simulation
%   % simulation can be started again from the last saved state of the 
%   % workspacem as set in "options.mcmc.i_save_workspace=10000;"
%   % Go to the folder with the mat-file and run "sippi_metropolis"
%   % with no input arguments, to restart sampling
%
% See also sippi_metropolis
%
%

% USE LESS DATA

% USE HIGHER UNCERTAINTY

% USE ADAPTIVE NORMALIZATION

%%

start_from_mat_file=0;

%% Figure out whether to continue sampling, or start crom scratch
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


if start_from_mat_file==0;
    % Start sampling from scratch
    mcmc.null='';
    if isfield(options,'mcmc');mcmc=options.mcmc;end
    
    %% NAME
    options.null='';
    if ~isfield(options,'txt');options.txt='';end
    if length(options.txt)>0
        options.txt=sprintf('%s_sippi_rejection_%s',datestr(now,'YYYYmmdd_HHMM'),options.txt);
    else
        options.txt=sprintf('%s_sippi_rejection',datestr(now,'YYYYmmdd_HHMM'));
    end
    
    
    %% INITIALIZE ASC FILE
    start_dir=pwd;
    try;
        mkdir(options.txt);
        cd(options.txt);
        addpath(['..',filesep])
        
        % copy training image file if used
        for im=1:length(prior);
            if isfield(prior{im},'ti')
                if ischar(prior{im}.ti)
                    try
                        if isunix
                            system(sprintf('cp ..%s%s . ',filesep,prior{im}.ti));
                        else
                            cmd=sprintf('copy ..%s%s ',filesep,prior{im}.ti);
                            system(cmd);
                        end
                    end
                end
            end
        end
    end
    for im=1:length(prior)
        filename_asc{im}=sprintf('%s_m%d%s',options.txt,im,'.asc');
        sippi_verbose(filename_asc{im},2);
        fid=fopen(filename_asc{im},'w');
        fclose(fid);
    end
    filename_mat=[options.txt,'.mat'];
    
    if ~isfield(mcmc,'T');mcmc.T=1;end
    
    if ~isfield(mcmc,'i_plot');mcmc.i_plot=500;end
    if ~isfield(mcmc,'adaptive_rejection');
        mcmc.adaptive_rejection=0;
    end
    if mcmc.adaptive_rejection==1
        if ~isfield(mcmc,'logLmax')&~isfield(mcmc,'Lmax')
            mcmc.logLmax=-inf;1e+300;
        end
    end
    if ~isfield(mcmc,'nite');mcmc.nite=1000;end
    if ~isfield(mcmc,'logLmax');
        if mcmc.adaptive_rejection==1;
            mcmc.logLmax=-100;
        else
            mcmc.logLmax=1;
        end
    end
    if isfield(mcmc,'Lmax');mcmc.Lmax=exp(mcmc.Lmax);end
    
    if ~isfield(mcmc,'i_save_workspace')
        mcmc.i_save_workspace=mcmc.i_plot*10;
    end
    
    
    prior=sippi_prior_init(prior);
    iacc=0;
    t0=now;
    
    mcmc.logL=zeros(1,mcmc.nite);
    
    if isfield(mcmc,'max_run_time_hours');
        mcmc.time_end = now + mcmc.max_run_time_hours/24;
    else
        mcmc.time_end = Inf;
    end
    
    sippi_verbose(sprintf('%s: STARTING rejection sampler in %s',mfilename,options.txt),-2)
    if mcmc.logLmax~=1;
        sippi_verbose(sprintf('%s: logLmax=%g',mfilename,mcmc.logLmax),-2)
    end
    if mcmc.adaptive_rejection==1;
        sippi_verbose(sprintf('%s: Adaptive rejection sampling, using logLmax=%g',mfilename,mcmc.logLmax),-2)
    end
    
end

i=0;
while i<=mcmc.nite;
    i=i+1;
    
    % OPTIONALLY LOAD FROM MAT_FILE
    if (i==1)&&(start_from_mat_file==1);
        % LOAD STATE FROM MATLAB
        load(mat_file);
    end
    
    
    % propose new model
    m_propose = sippi_prior(prior);
    [d,forward,prior]=sippi_forward(m_propose,forward,prior);
    [logL,L,data]=sippi_likelihood(d,data);
    
    logLPacc = (1./mcmc.T).*(logL-mcmc.logLmax);
    
    if log(rand(1))<logLPacc
        sippi_verbose(sprintf('%s: %06d/%06d ACCEPT logLPacc=%4.1g, Pacc=%4.1g',mfilename,i,mcmc.nite,logLPacc,exp(logLPacc)),1);
        iacc=iacc+1;
        mcmc.logL(iacc)=logL;
        for im=1:length(prior)
            fid=fopen(filename_asc{im},'a+');
            fprintf(fid,' %10.7g ',m_propose{im}(:));
            fprintf(fid,'\n');
            fclose(fid);
        end
        
    end
    if (i/mcmc.i_plot)==round(i/mcmc.i_plot)
        
        [t_end_txt,t_left_seconds]=time_loop_end(t0,i,mcmc.nite);
        nite=mcmc.nite;
        
        % time left if using
        if ~isinf(mcmc.time_end)
            t_left_seconds_time_end = 3600*24*(mcmc.time_end-now);
            if (t_left_seconds_time_end<t_left_seconds)
                t_end_txt = datestr(mcmc.time_end);%'time_limit';
                % compute reamining number of iterations
                time_per_ite = ((now-t0)/i);
                i_left = (mcmc.time_end-now)/time_per_ite;
                nite = i+ceil(i_left);
            end
        end
        sippi_verbose(sprintf('%s: %06d/%06d (%10s) nacc=%06d - %s',mfilename,i,nite,t_end_txt,iacc),-1)
        
    end
    
    %% ADAPTIVE REJECTION
    if (mcmc.adaptive_rejection==0)
        % Traditional rejection sampling
        
    else
        % Adaptive rejection sampling
        if logL>mcmc.logLmax
            sippi_verbose(sprintf('%s: i=%06d,  new log(maxL) = %g (%g)',mfilename,i,logL,mcmc.logLmax))
            mcmc.logLmax=logL;
        end
        
    end
    
    %% CHECK FOR TIME LIMIT
    if (now>mcmc.time_end);
        sippi_verbose(sprintf('%s: i=%06d,  TIME LIMIT REACHED!',mfilename,i))
        
        break
        %else
        %    disp(sprintf(' %s - %s',datestr(now),datestr(mcmc.time_end)))
    end
    
    %% SAVE WORKSPACE
    if ((i/(mcmc.i_save_workspace))==round( i/(mcmc.i_save_workspace) ))
        try
            sippi_verbose(sprintf('%s: i=%d, saving workspace to ''%s''',mfilename,i,filename_mat))
            save(filename_mat,'-v7.3')
        catch
            sippi_verbose(sprintf('%s: failed to save data to %s',mfilename,filename_mat))
        end
    end
    
    
    
end
mcmc.logL=mcmc.logL(1:iacc);

options.mcmc=mcmc;

save(filename_mat)
sippi_verbose(sprintf('%s : DONE rejection sampling in %s',mfilename,options.txt),-2)


%%
cd(start_dir);




end
