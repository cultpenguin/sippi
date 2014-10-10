% sippi_prior_dsim : Direct simulation in SIPPI
%
%
% TMH/2014
%
% See also: sippi_prior_init, sippi_prior
%
function [m_propose,prior]=sippi_prior_dsim(prior,m_current,ip);

if ~exist('dsim.m','file')
    disp(sptinf('%s: dsim.m is not in the path!!',mfilename));
    return
end

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if ~isfield(prior{ip},'ti')
    d=load('topography');
    ti=flipud(d.topo);
    ti(find(ti>0))=1;
    ti(find(ti<=0))=0;
    prior{ip}.ti=ti;
end



prior{ip}.options.null='';

if nargin == 1
    if prior{ip}.ndim==1;
        SIM_data=NaN.*ones(1,length(prior{ip}.x));
    elseif prior{ip}.ndim==2;
        SIM_data=NaN.*ones(length(prior{ip}.x),length(prior{ip}.x));
    else
        SIM_data=NaN;
        disp(sprintf('%s: 3D not supported for DSIM yet',mfilename))
    end
else
    % PERTURB
    SIM_data=m_current{ip};
    N=prod(size(SIM_data));
    i_resample=randomsample(N,prior{1}.seq_gibbs.step*N);
   SIM_data(i_resample)=NaN;
end
[m_propose{ip},prior{ip}.options]=dsim(prior{ip}.ti,SIM_data,prior{ip}.options);

%m_propose{ip}=normcdf(prior{ip}.randn,0,1)*(prior{ip}.max-prior{ip}.min)+prior{ip}.min;
    
