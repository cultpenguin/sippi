% sippi_prior_dummy : 'Dummy' prior for SIPPI. 
%
%% Example:
%   ip=1;
%   prior{ip}.type='dummy';
%   sippi_prior(prior)
%       and
%
function [m_propose,prior]=sippi_prior_dummy(prior,m_current,ip);

if nargin<3;
    ip=1;
end

m_propose{ip}=NaN;