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
%   options.mcmc.logLmax
%
%   options.mcmc.rejection_normalize_log = log(options.mcmc.Lmax)
%
%   options.mcmc.adaptive_rejection=1, adaptive setting of maxiumum likelihood
%                  (def=[0])
%                  At each iteration Lmax will be set if log(L(m_cur)=>options.mcmc.logLmax
%
%
%   options.mcmc.max_run_time_hours = 1; % maximum runtime in hours
%                                        % (overrides options.mcmc.nite if needed)
%
% See also sippi_metropolis
%
%

% USE LESS DATA

% USE HIGHER UNCERTAINTY

% USE ADAPTIVE NORMALIZATION


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
for im=1:length(prior)
    filename_asc{im}=sprintf('%s_m%d%s',options.txt,im,'.asc');
    disp(filename_asc{im});
    fid=fopen(filename_asc{im},'w');
    fclose(fid);
end
filename_mat=[options.txt,'.mat'];

mcmc.null='';
if isfield(options,'mcmc');mcmc=options.mcmc;end


if ~isfield(mcmc,'i_plot');mcmc.i_plot=500;end
if ~isfield(mcmc,'adaptive_rejection');mcmc.adaptive_rejection=0;end
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
if ~isfield(mcmc,'rejection_normalize_log');mcmc.rejection_normalize_log = mcmc.logLmax;end

prior=sippi_prior_init(prior);
iacc=0;
t0=now;

mcmc.logL=zeros(1,mcmc.nite);

if isfield(mcmc,'max_run_time_hours');
    mcmc.time_end = now + mcmc.max_run_time_hours/24;
else
    mcmc.time_end = Inf;
end

for i=1:mcmc.nite
    
    m_propose = sippi_prior(prior);
    if isfield(forward,'forward_function');
        [d,forward,prior,data]=feval(forward.forward_function,m_propose,forward,prior,data);
    else
        [d,forward,prior,data]=sippi_forward(m_propose,forward,prior,data);
    end
    [logL,L,data]=sippi_likelihood(d,data);
    
    logLPacc = logL-mcmc.rejection_normalize_log;
    
    if log(rand(1))<logLPacc
        iacc=iacc+1;
        mcmc.logL(iacc)=logL;
        %disp(sprintf('%s : i=%04d, nacc=%04d',mfilename,i,iacc));
        % write current model to disc
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
        disp(sprintf('%s : %06d/%06d (%10s) nacc=%06d - %s',mfilename,i,nite,t_end_txt,iacc))
        
    end
    if ((i/(10*mcmc.i_plot))==round( i/(10*mcmc.i_plot) ))
        save(filename_mat)
    end
    
    if (mcmc.adaptive_rejection==0)
        % Traditional rejection sampling
        
    else
        % Adaptive rejection sampling
        if logL>mcmc.rejection_normalize_log
            disp(sprintf('%s : i=%06d,  new log(maxL) = %g (%g)',mfilename,i,logL,mcmc.rejection_normalize_log))
            mcmc.rejection_normalize_log=logL;
        end
        
    end
    
    if (now>mcmc.time_end);
        break
        %else
        %    disp(sprintf(' %s - %s',datestr(now),datestr(mcmc.time_end)))
    end
    
    
end
mcmc.logL=mcmc.logL(1:iacc);

options.mcmc=mcmc;

save(filename_mat)
disp(sprintf('%s : DONE rejection on %s',mfilename,options.txt))


%%
cd(start_dir);




end
