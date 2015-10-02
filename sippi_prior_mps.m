% sippi_prior_mps : prior based on MPS
%
%                      Using SNESIM/ENESIM FROM 
%                      https://github.com/cultpenguin/mps
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
function [m_propose,prior]=sippi_prior_mps(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if ~isfield(prior{ip},'ti')
    prior{ip}.ti=channels;
end
%prior{ip}.S.method='mps_snesim_list';
prior{ip}.S.method='mps_snesim_tree';
%prior{ip}.S.template_size=[9 9 1];
prior{ip}.S.nreal=1;
%prior{ip}.S.n_multiple_grids=3;
%prior{ip}.S.shuffle_simulation_grid=1;
prior{ip}.S.parameter_filename='uncond.txt';

if prior{ip}.ndim==1;
    SIM=zeros(prior{ip}.dim(1));    
elseif prior{ip}.ndim==1;
    SIM=zeros(prior{ip}.dim(2),prior{ip}.dim(1));    
else
    SIM=zeros(prior{ip}.dim(2),prior{ip}.dim(1),prior{ip}.dim(3));    
end

% RUN FORWARD
[m_propose{ip},prior{ip}.S]=mps_cpp(prior{ip}.ti,SIM,prior{ip}.S);



