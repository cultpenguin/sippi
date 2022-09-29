% sippi_prior_gamma : gamma prior for SIPPI
%
% See also: sippi_prior_init, sippi_prior
%
function [m_propose,prior]=sippi_prior_gamma(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if ~isfield(prior{ip},'rprior')
    % This can be any prior type that returns values between 0 and 1!
    prior{ip}.rprior{1}.type='uniform';
    prior{ip}.rprior{1}.min=0;
    prior{ip}.rprior{1}.max=1;
    prior{ip}.rprior{1}.x=prior{ip}.x;
    prior{ip}.rprior{1}.y=prior{ip}.y;
    prior{ip}.rprior{1}.z=prior{ip}.z; 
end

if ~isfield(prior{ip},'a');prior{ip}.a=10;end
if ~isfield(prior{ip},'b');prior{ip}.b=2;end


if nargin==1;
    [prior{ip}.rand,prior{ip}.rprior]=sippi_prior(prior{ip}.rprior);
else
    % update step
    prior{1}.rprior{1}.seq_gibbs=prior{1}.seq_gibbs;
    [prior{ip}.rand,prior{ip}.rprior]=sippi_prior(prior{ip}.rprior,prior{ip}.rand);
end
m_propose{ip}=gaminv(prior{ip}.rand{1},prior{ip}.a,prior{ip}.b);
    

