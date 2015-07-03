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
%    % optionally a specific random seed can be set using
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
%
% ip=1;
% prior{ip}.type='visim';
% prior{ip}.d_target=d_target;
% prior{ip}.cax=[min(d_target) max(d_target)];
% prior{ip}.x=1:1:80;
% prior{ip}.y=1:1:80;
% prior{ip}.Cm=sprintf('%g Gau(20)',var(d_target));
% [m,prior]=sippi_prior(prior);
% sippi_plot_prior(prior,m);
%
%
% See also: sippi_prior, visim
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

%% RANDOM SEED
prior{ip}.V.cond_sim=0;
if isfield(prior{ip},'seed');
    prior{ip}.V.rseed=prior{ip}.seed;
    elsec
    prior{ip}.V.rseed=ceil(rand(1).*1e+6);
    sippi_verbose(sprintf('%s : setting seed (%d)for VISIM',mfilename,prior{ip}.V.rseed),2)
end

%% SET TARGET DISTRIBUTION
if isfield(prior{ip},'d_target')
    f_cond=sprintf('d_target_%02d.eas',ip);
    if ~exist(f_cond,'file');
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
    
    prior{ip}.V.gmean=mean(prior{ip}.d_target);
    
end

%% SEQUENTIAL GIBBS
if nargin>1
    % SEQUENTIAL GIBBS
    mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
    %prior{ip}.S=sgems_set_resim_data(prior{ip}.S,m_current,prior{ip}.seq_gibbs.step,prior{ip}.seq_gibbs.type);
    [prior{ip}.V, i_resim]=visim_set_resim_data(prior{ip}.V,m_current{ip},prior{ip}.seq_gibbs.step,[],[],prior{ip}.seq_gibbs.type);
    prior{ip}.seq_gibbs.used=i_resim;
    if isempty(i_resim)
        prior{ip}.V.cond_sim=0;
    else
        prior{ip}.V.cond_sim=2;
    end
    
end

%% RUN VISIM
prior{ip}.V=visim(prior{ip}.V);
m_propose{ip} = prior{ip}.V.D';



