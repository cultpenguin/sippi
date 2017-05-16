% sippi_prior_mixsim : MIXsim prior for SIPPI (only 2D)
%
% Example:  
%  ip=1; 
%  prior{ip}.type='mixsim'; % MIXSIM type
%  prior{ip}.x=[1:1:20]; % X array 
%  prior{ip}.y=[1:1:40]; % Y array 
%  TI=channels;
%  TI=TI(3:3:100,1:1:10)+1;
%  prior{1}.TI=TI; % must be integer values starting with 1! (no zero values)
%  m=sippi_prior(prior);
%
%  sippi_plot_prior_sample(prior);
%
% See also sippi_prior, mixsim_2D
% 

%
% Please refer to the following paper when using this algorithm:
% Cordua, K.S., T.M. Hansen, M.L. Gulbrandsen, C. Barnes, and
% K. Mosegaard, 2016. Mixed-point geostatistical simulation: A
% combination of two- and multiple-point geostatistics. Geophysical
% Research Letters.
%
% K.S. Cordua and T.M. Hansen (2016)
%
% See also: sippi_sequential_gibbs_resim
%
function [m_propose,prior]=sippi_prior_mixsim(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if ~isfield(prior{ip},'TI');
    sippi_verbose(sprintf('%s: No TI field (training image) set!',mfilename))
    m_propose=[];
    return
end


prior{ip}.options.do_cond=0; % no conditional simulation
%% Sequential gibbs resampling
if nargin>1
    prior{ip}.options.do_cond=1;
    d_cond=sippi_get_resim_data(m_current,prior,ip);
    % set hard data
    if ~isempty(d_cond);
        prior{ip}.options.data=[d_cond(:,4), d_cond(:,2), d_cond(:,1)];
    end
end


%%
prior{ip}.options.null=[];
[m_propose{ip},prior{ip}.options]=mixsim_2D(prior{ip}.dim(2),prior{ip}.dim(1),1,prior{ip}.TI,prior{ip}.options);

debug=0;
if (debug==1)
    if nargin>1
        figure(1);clf;
        subplot(1,3,1);imagesc(m_current{ip});axis image;
        try
            subplot(1,3,2);scatter(prior{ip}.options.data(:,3),prior{ip}.options.data(:,2),20,prior{ip}.options.data(:,1),'filled');axis image
            set(gca,'ydir','revers')
        end
        subplot(1,3,3);imagesc(m_propose{ip});axis image;
        drawnow;
    end
end