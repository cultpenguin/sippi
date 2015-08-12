% sippi_prior_snesim : SNESIM type Gaussian prior for SIPPI
%
%                      Using SNESIM form 
%                      https://github.com/SCRFpublic/snesim-standalone
% 
%% Example:
%    ip=1;
%    prior{ip}.type='snesim';
%    prior{ip}.x=1:1:80;
%    prior{ip}.y=1:1:80;
%    prior{ip}.ti=channels;
%    % prior{ip}.ti=maze;
%
%    m=sippi_prior(prior);
%    sippi_plot_prior(prior,m)
%    figure(1);imagesc(prior{ip}.ti);axis image
%
%% Example: scaling and rotation
%    ip=1;
%    prior{ip}.type='snesim';
%    prior{ip}.x=1:1:80;
%    prior{ip}.y=1:1:80;
%    prior{ip}.ti=channels;
%    prior{ip}.scaling=[.1];
%    prior{ip}.rotation=[10];
%
%    m=sippi_prior(prior);
%    sippi_plot_prior(prior,m)
%    figure(1);imagesc(prior{ip}.ti);axis image
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
% See also: sippi_prior, ti
%
function [m_propose,prior]=sippi_prior_snesim(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end


% SNESIM

%% REMOVE CONDITIONAL DATA.
%% FIX : NEED TO CHANGE TO HANDLE CONDITIONAL DATA
%if isfield(prior{ip}.S,'f_obs')
%    prior{ip}.S=rmfield(prior{ip}.S,'f_obs');
%end
%prior{ip}.S.XML.parameters.Hard_Data.grid='';
%prior{ip}.S.XML.parameters.Hard_Data.property='';

% force nsim=1 in sippi
prior{ip}.S.nsim=1;

% set random seed
if isfield(prior{ip},'seed');
    prior{ip}.S.rseed.seed;
else
    prior{ip}.S.rseed=ceil(rand(1).*1e+6);
end

% optionally set rotation and affinity
set_aff=0;
if isfield(prior{ip},'rotation')|isfield(prior{ip},'scaling')
  set_aff=1;
end
if set_aff==1
  if ~isfield(prior{ip},'rotation'), prior{ip}.rotation=1; end
  if ~isfield(prior{ip},'scaling'), prior{ip}.scaling=1; end
  prior{ip}.S=snesim_set_rotation_affinity(prior{ip}.S,prior{ip}.rotation,1./prior{ip}.scaling);
  prior{ip}.S.frotaff.use=1;
end

if nargin>1
    % SEQUENTIAL GIBBS
    if isfield(prior{ip},'index_values');
        m = zeros(size(m_current{ip}))-1;
        for i=1:length(prior{ip}.index_values)
            try
                m(find(m_current{ip}==prior{ip}.m_values(i)))=prior{ip}.index_values(i);
            end
        end
        m_current{ip}=m;
    end
    sippi_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
    %prior{ip}.S=snesim_set_resim_data(prior{ip}.S,prior{ip}.S.D,[10 10]);
    
    prior{ip}.S=snesim_set_resim_data(prior{ip}.S,m_current{ip},prior{ip}.seq_gibbs.step,prior{ip}.seq_gibbs.type);

end

%% RUN SNESIM
prior{ip}.S = snesim(prior{ip}.S,prior{ip}.x,prior{ip}.y,prior{ip}.z);
m_propose{ip} = prior{ip}.S.D;

if isfield(prior{ip},'index_values');
    for i=1:length(prior{ip}.index_values)
        m_propose{ip}(find(m_propose{ip}==prior{ip}.index_values(i)))=prior{ip}.m_values(i);
    end
end
