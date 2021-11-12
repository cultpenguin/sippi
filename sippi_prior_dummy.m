% sippi_prior_dummy : 'Dummy' prior for SIPPI. 
%
%% Example:
%   ip=1;
%   prior{ip}.type='dummy';
%   sippi_prior(prior)
%
function [m_propose,prior]=sippi_prior_dummy(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if nargin>1
    m_propose{ip}=m_current{ip};
else
    m_propose{ip}=NaN;
end