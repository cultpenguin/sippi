% sippi_prior_visim : VISIM type Gaussian prior for SIPPI
%
%% Example:
%    ip=1;
%    prior{ip}.type='visim';
%    prior{ip}.x=1:1:80;
%    prior{ip}.y=1:1:80;
%    prior{ip}.Cm='1 Sph(60)';
%    m=sippi_prior(prior);
%    sippi_plot_prior(prior,m)
%
%    % optionally a specific random can be set using
%    prior{ip}.seed=1;
%
%% Sequential Gibbs sampling type 1 (box selection of pixels)
%    prior{ip}.seq_gibbs.type=1;
%    prior{ip}.seq_gibbs.step=10; % resim data in 10x10 pixel grids
%    prior{ip}.cax=[-2 2];
%    [m,prior]=sippi_prior(prior);
%    for i=1:1000;
%       [m,prior]=sippi_prior(prior,m);
%       sippi_plot_prior(prior,m);
%       drawnow;
%    end
%
%% Sequential Gibbs sampling type 2 (random pixels)
%    prior{ip}.seq_gibbs.type=2;%
%    prior{ip}.seq_gibbs.step=.6; % Resim 60% of data
%    [m,prior]=sippi_prior(prior);
%    for i=1:10;
%       [m,prior]=sippi_prior(prior,m);
%       sippi_plot_prior(prior,m);
%       drawnow;
%    end
%
%% TARGET DISTRIBUTION
% clear prior
% d_target=[7 8 9 10 11 11 12];
% d_target=[7 8 9 10 14 15 20];
%
% ip=1;
% prior{ip}.type='visim';
% prior{ip}.method='sgsim';
% % prior{ip}.method='dssim';
% prior{ip}.d_target=d_target;
% prior{ip}.cax=[min(d_target) max(d_target)];
% prior{ip}.x=1:1:80;
% prior{ip}.y=1:1:80;
% prior{ip}.Cm=sprintf('%g Gau(20)',var(d_target));
% prior{ip}.Cm=sprintf('%g Gau(20)',1);
% [m,prior]=sippi_prior(prior);
% sippi_plot_prior(prior,m);
%
% prior{ip}.seq_gibbs.step=16;
% prior{ip}.seq_gibbs.type=1;
% for i=1:10;
%    [m,prior]=sippi_prior(prior,m);
%    sippi_plot_prior(prior,m);
%    drawnow
% end
%
%
% See also: sippi_prior, visim, nscore, inscore
%
function [m_propose,prior]=sippi_prior_visim(prior,m_current,ip);
if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

% VISIM PRIOR
if isfield(prior{ip},'Va');
    % update VISIM covariance settings
    if ~isstruct(prior{ip}.Va);
        Va=deformat_variogram(prior{ip}.Va);
    else
        Va=prior{ip}.Va;
    end
    [prior{ip}.V]=visim_set_variogram(prior{ip}.V,prior{ip}.Va);
end

if isfield(prior{ip},'m0')
    prior{ip}.V.gmean=prior{ip}.m0; % set global mean
end

% SIMULATION METHOD
if ~isfield(prior{ip},'method');
    prior{ip}.method='dssim';
%    if isfield(prior{ip},'d_target')
%        %prior{ip}.method='sgsim';
%    else
%        prior{ip}.method='sgsim';
%    end
end

%% RANDOM SEED
prior{ip}.V.cond_sim=0;
if isfield(prior{ip},'seed');
    prior{ip}.V.rseed=prior{ip}.seed;
else
    prior{ip}.V.rseed=ceil(rand(1).*1e+6);
    sippi_verbose(sprintf('%s : setting seed (%d)for VISIM',mfilename,prior{ip}.V.rseed),2)
end

%% CONDITIONAL POINT DATA, d_obs
if isfield(prior{ip},'d_obs')
    
    if size(prior{ip}.d_obs,2)==4
        prior{ip}.d_obs(:,5)=0;
    end
    
    i_hard=find(prior{ip}.d_obs(:,5)==0);
    i_soft=find(prior{ip}.d_obs(:,5)~=0);
    
    sippi_verbose(sprintf('%s: Using %d hard data',mfilename,length(i_hard)),-10);
    sippi_verbose(sprintf('%s: Using %d soft data',mfilename,length(i_soft)),-10);
    
    useHardPoint=0;
    useSoftPoint=0;
    
    prior{ip}.V.fconddata.fname='d_obs.eas';
    if isempty(i_hard)
        delete(prior{ip}.V.fconddata.fname)
    else
        useHardPoint=1;        
        write_eas(prior{ip}.V.fconddata.fname,prior{ip}.d_obs(i_hard,1:4));
        prior{ip}.V.cond_sim=2; % only point data
    
    end
    if ~isempty(i_soft)
        useSoftPoint=1;
        % volume data
        
        prior{ip}.V.fvolgeom.fname='d_volgeom.eas';
        prior{ip}.V.fvolsum.fname='d_volsum.eas';
        clear d_volgeom d_volsum
        for i=1:length(i_soft)
            d_volgeom(i,:)=[prior{ip}.d_obs(i_soft(i),1:3) i 1];
            d_volsum(i,:)= [i 1 prior{ip}.d_obs(i_soft(i),4:5)];
        end
        write_eas(prior{ip}.V.fvolgeom.fname,d_volgeom);
        write_eas(prior{ip}.V.fvolsum.fname,d_volsum);
        if useHardPoint==1
            prior{ip}.V.cond_sim=1; % hard and soft data
        else
            prior{ip}.V.cond_sim=3; % soft data
        end
        
        
        prior{ip}.V.volnh.method=2; %--> BAD BAD RESULTS
        prior{ip}.V.densitypr=0;
        prior{ip}.V.debuglevel=-1;
        
    end
end


%% SET TARGET DISTRIBUTION
if isfield(prior{ip},'d_target')
    if strcmp(prior{ip}.method,'dssim');
        
        %% DSSIM
        % make sure each visim type has different filename for target dist.       
        if isfield(prior{ip},'visim_id');
            f_cond=sprintf('d_target_%02d.eas',prior{ip}.visim_id);
        else
            f_cond=sprintf('d_target.eas',ip);
        end
        if ~exist([pwd,filesep,f_cond],'file');
            write_eas(f_cond,prior{ip}.d_target(:));
            sippi_verbose(sprintf('%s : writing target distribution to %s',mfilename,f_cond));
        end
        % set data the define a realization of a 1D marginal distribution
        prior{ip}.V.refhist.fname=f_cond;
        % Use the 1D marginal distribution as terget distritbution
        prior{ip}.V.ccdf=1;
        % Treat the 1D marginal distribution as a discrete distribution
        % in this case only actual values from the 1D marginal will
        % be realized
        prior{ip}.V.refhist.do_discrete=1;
        
        
        if isfield(prior{ip},'min')
            prior{ip}.V.tail.zmin=prior{ip}.min;
        else
            prior{ip}.V.tail.zmin=min(prior{ip}.d_target);
        end
        if isfield(prior{ip},'max')
            prior{ip}.V.tail.zmin=prior{ip}.max;
        else
            prior{ip}.V.tail.zmax=max(prior{ip}.d_target);
        end
        
        % update global mean and variance in visim parameter file!
        sippi_verbose(sprintf('%s: Updating global mean and variance from target dist',mfilename),10);
        
        prior{ip}.V.gmean=mean(prior{ip}.d_target);
        if ~isstruct(prior{ip}.Va);
            prior{ip}.Va=deformat_variogram(prior{ip}.Va);
        end
        Va_par=prior{ip}.Va;
        gvar_Va=sum([Va_par.par1]);
        gvar_d_target=var(prior{ip}.d_target);
        prior{ip}.V.gvar=gvar_d_target;
        
        for j=1:length(prior{ip}.Va);
            prior{ip}.Va(j).par1 = prior{ip}.Va(j).par1 * (gvar_d_target./gvar_Va);
        end
        [prior{ip}.V]=visim_set_variogram(prior{ip}.V,prior{ip}.Va);
       
        prior{ip}.m0= prior{ip}.V.gmean;
        
        
    else
        %% SGSIM
        % setup normal score transform
        %if (isfield(prior{ip},'d_target'))&(~isfield(prior{ip},'o_nscore'))
        if (~isfield(prior{ip},'o_nscore'))&&(isfield(prior{ip},'d_target'));
            % UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION
            d_min=min(prior{ip}.d_target);
            d_max=max(prior{ip}.d_target);
            [d_nscore,o_nscore]=nscore(prior{ip}.d_target,1,1,d_min,d_max,0);
            prior{ip}.o_nscore=o_nscore;
            
            % force mean to zero and variance to 1
            sippi_verbose(sprintf('%s: Updating global mean and variance to N(0,1)',mfilename),10);
            if ~isstruct(prior{ip}.Va);
                prior{ip}.Va=deformat_variogram(prior{ip}.Va);
            end
            Va_par=prior{ip}.Va;
            gvar=sum([Va_par.par1]);
 
            for j=1:length(prior{ip}.Va);
                prior{ip}.Va(j).par1 = prior{ip}.Va(j).par1./gvar;
            end
            prior{ip}.V.ccdf=0; % DO NOT USE TARGET DISTRIBUTION
            prior{ip}.V.gvar=1;
            prior{ip}.V.gmean=0;
            %prior{ip}.V.tail.zmin=-5;
            %prior{ip}.V.tail.zmin=5;
            prior{ip}.m0=0;
            
        
        end
     
    end
    
end

%% SEQUENTIAL GIBBS
if nargin>1
    
    if ~isfield(prior{ip}.seq_gibbs,'pos')
        prior{ip}.seq_gibbs.pos=[];
    end
    
    % SEQUENTIAL GIBBS
    mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
    
    m=m_current{ip};
    
    % if using SGSIM methods and using target distribution, then perform
    % forward normal score! Thius will add some uncertainty
    if (strcmp(prior{ip}.method,'sgsim'))&(isfield(prior{ip},'o_nscore'))
        
        %if isfield(prior{ip}.V,'D');
        %    m=prior{ip}.V.D';
        %else
            m=nscore(m,prior{ip}.o_nscore);
        %end
    end
    
    %[prior{ip}.V, i_resim]=visim_set_resim_data(prior{ip}.V,m,prior{ip}.seq_gibbs.step,[],[],prior{ip}.seq_gibbs.type);
    [prior{ip}.V, i_resim]=visim_set_resim_data(prior{ip}.V,m,prior{ip}.seq_gibbs.step,prior{ip}.seq_gibbs.pos,[],prior{ip}.seq_gibbs.type);
    prior{ip}.seq_gibbs.used=i_resim;
    if isempty(i_resim)
        prior{ip}.V.cond_sim=0;
    else
        prior{ip}.V.cond_sim=2;
    end
    
end

%% RUN VISIM
[prior{ip}.V,status,result]=visim(prior{ip}.V);
% Check that output has been generated
if ~isfield(prior{ip}.V,'D')
    disp(sprintf('%s: Something went running visim on "%s"',mfilename,'visim.par'));
    disp(result)
    m_propose{ip}=[];
    return
end

m_propose{ip} = prior{ip}.V.D';

%% PERFORM NORMAL SCORE OF NEEDED
if (strcmp(prior{1}.method,'sgsim'))&(isfield(prior{ip},'o_nscore'))
    if ~isstruct(prior{ip}.Va);
        prior{ip}.Va=deformat_variogram(prior{ip}.Va);
    end
    Va_par=prior{ip}.Va;
    gvar=sum([Va_par.par1]);
    m_propose{ip}=m_propose{ip}./sqrt(gvar);
    m_propose{ip}=inscore(m_propose{ip},prior{ip}.o_nscore);
    
    
    % add mean model
    m_propose{ip}=m_propose{ip}+prior{ip}.m0;
    
end
     



