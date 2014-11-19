% sippi_prior_sisim: SISIM (SGeMS) type prior for SIPPI
%
%% Example:
%    ip=1;
%    prior{ip}.type='sisim';
%    prior{ip}.x=1:1:80;
%    prior{ip}.y=1:1:80;
%    prior{ip}.Cm='1 Sph(60)';
%    prior{ip}.marginal_prob=[.1 .4 .5];
%    m=sippi_prior(prior);
%    sippi_plot_prior(prior,m)
%
%    % optionally a specific random seed can be set using
%    prior{ip}.seed=1;
%
%% Sequential Gibbs sampling type 1 (box selection of pixels)
%    prior{ip}.seq_gibbs.type=1;%
%    prior{ip}.seq_gibbs.step=10; % resim data in 10x10 pixel grids
%    [m,prior]=sippi_prior(prior);
%    for i=1:10;
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
% See also: sippi_prior, sgems
%
function [m_propose,prior]=sippi_prior_sisim(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if ~isfield(prior{ip},'S');
    prior{ip}.S=sgems_get_par('sisim');
end
if isfield(prior{ip},'marginal_prob');
    prior{ip}.S.XML.parameters.Marginal_Probabilities.value=prior{ip}.marginal_prob;
    % prior{ip}.marginal_prob=ones(1,n_ind)/n_ind;
    n_ind=length(prior{ip}.S.XML.parameters.Marginal_Probabilities.value);
    prior{ip}.S.XML.parameters.Nb_Indicators.value=n_ind;
end


if isfield(prior{ip},'Va');
    va_xml=sgems_variogram_xml(prior{ip}.Va);
    prior{ip}.S.XML.parameters.Variogram_Median_Ik=va_xml;
end

% NREAL = 1
prior{ip}.S.XML.parameters.Nb_Realizations.value=1;
% SEED
if isfield(prior{ip},'seed');
    prior{ip}.S.XML.parameters.Seed.value=prior{ip}.seed;
else
    prior{ip}.S.XML.parameters.Seed.value=ceil(rand(1).*1e+6);
end

% REMOVE CONDITIONAL DATA.
% FIX : NEED TO CHANGE TO HANDLE CONDITIONAL DATA
if isfield(prior{ip}.S,'f_obs')
    prior{ip}.S=rmfield(prior{ip}.S,'f_obs');
end
prior{ip}.S.XML.parameters.Hard_Data_Grid.value='';
prior{ip}.S.XML.parameters.Hard_Data_Property.value='';

% SEQ GIBBS
if nargin>1
    mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
    
    prior{ip}.S=sgems_set_resim_data(prior{ip}.S,m_current{ip},prior{ip}.seq_gibbs.step,prior{ip}.seq_gibbs.type);
    %[data,used]=sgems_set_resim_data(prior{ip}.S,m_current{ip},prior{ip}.seq_gibbs.step,prior{ip}.seq_gibbs.type);
    %prior{ip}.S=data;
    %prior{ip}.seq_gibbs.used=used;
end

% SISIM SIM
prior{ip}.S = sgems_grid(prior{ip}.S);
m_propose{ip} = prior{ip}.S.D';
