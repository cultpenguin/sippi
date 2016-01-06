% sippi_prior_plurigaussian: Plurigaussian type prior for SIPPI
%
%% Example:
%   % PluriGaussian based on one Gaussian model / truncated Gaussian
%   ip=1;
%   prior{ip}.type='plurigaussian';
%   prior{ip}.x=1:1:80;
%   prior{ip}.y=1:1:80;
%   prior{ip}.Cm='.001 Nug(0) + 1 Gau(10)';
%   prior{ip}.pmap=[0 0 0 0 1 1 0 0 2 2 2];
%   % PluriGaussian based on two Gaussian models 
%   ip=ip+1;
%   prior{ip}.type='plurigaussian';
%   prior{ip}.x=1:1:80;
%   prior{ip}.y=1:1:80;
%   prior{ip}.Cm{1}='1 Nug(0) + 1 Gau(10)';
%   prior{ip}.Cm{2}='1 Nug(0) + 1 Gau(10)';
%   prior{ip}.pmap=[0 0 0 1 1; 1 2 0 0 1 1; 1 1 1 1 1];
%
%   [m,prior]=sippi_prior(prior);
%   sippi_plot_prior(prior,m);
%
%% Sequential Gibbs sampling
%   prior{1}.seq_gibbs.step=.1;
%   for i=1:100;
%       [m,prior]=sippi_prior(prior,m);
%       sippi_plot_prior(prior,m);
%       caxis([8 12]);drawnow;
%   end
%
% See also: sippi_prior
%
function [m_propose,prior]=sippi_prior_plurigaussian(prior,m_current,ip)

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if (isfield(prior{ip},'d_target'))&(~isfield(prior{ip},'o_nscore'))
    % UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION
    d_min=min(prior{ip}.d_target);
    d_max=max(prior{ip}.d_target);
    [d_nscore,o_nscore]=nscore(prior{ip}.d_target,1,1,d_min,d_max,0);
    prior{ip}.o_nscore=o_nscore;
end

%% Remove Rand Number if they are not needed
if (nargin==1)&isfield(prior{ip},'z_rand');
    prior{ip}=rmfield(prior{ip},'z_rand');
end


%% CHECK FORMAT

% test if simple 1D prior is set 
if isfield(prior{ip},'Cm')
    % One Gaussian 
    if ~isfield(prior{ip},'pg_prior')
        prior{ip}.pg_prior{1}.type='fftma';
        prior{ip}.pg_prior{1}.x=prior{ip}.x;
        prior{ip}.pg_prior{1}.y=prior{ip}.y;
        prior{ip}.pg_prior{1}.z=prior{ip}.z;
        prior{ip}.pg_prior{1}.m0=1;
        prior{ip}.pg_prior{1}.Cm=prior{ip}.Cm;
    end
end
        
prior{ip}.NG=length(prior{ip}.pg_prior);

% TEST
for ipg=1:prior{ip}.NG
    if ~isfield(prior{ip}.pg_prior{ipg},'x')
        prior{ip}.pg_prior{ipg}.x=prior{ip}.x;
    end
    if ~isfield(prior{ip}.pg_prior{ipg},'y')
        prior{ip}.pg_prior{ipg}.y=prior{ip}.y;
    end
    if ~isfield(prior{ip}.pg_prior{ipg},'z')
        prior{ip}.pg_prior{ipg}.z=prior{ip}.z;
    end
    if ~isfield(prior{ip}.pg_prior{ipg},'m0')
        prior{ip}.pg_prior{ipg}.m0=0;
    end
    if ~isfield(prior{ip}.pg_prior{ipg},'Cm')
        prior{ip}.pg_prior{ipg}.Cm=prior{ip}.Cm;
    end
    if ~isfield(prior{ip}.pg_prior{ipg},'type')
        prior{ip}.pg_prior{ipg}.type='fftma';
    end
end

%% SIMULATE PLURIGAUSSIAN DATA
if nargin>1
    % Sequential Gibbs?
    % [pg_m_propose,prior{1}.pg_prior]=sippi_prior(prior{1}.pg_prior);
    disp('STOP')
    keyboard
else
    % unconditional simulation
    [pg_m_propose,prior{ip}.pg_prior]=sippi_prior(prior{ip}.pg_prior);
end

%% CONVERT N(0,1) using PG_MAP
%m_propose = plurigaussian_transform(pg_m_propose,prior{ip}.pg_map);
m_propose{ip} = pg_m_propose{1}.*0;
m_propose{ip}(find(pg_m_propose{1}>-1))=1;
m_propose{ip}(find(pg_m_propose{1}>0))=0;
m_propose{ip}(find( (pg_m_propose{1}>1.5) ))=3;
m_propose{ip}(find( (pg_m_propose{1}>1.5)&(pg_m_propose{2}>1.5) ))=4;


       