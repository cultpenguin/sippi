% sippi_prior_dsim : Direct simulation in SIPPI
%
% Example: 
%  
% prior{1}.type='dsim';
% prior{1}.x=1:1:40;;
% prior{1}.y=1:1:30;;
% prior{1}.ti=channels;;
%
% m=sippi_prior(prior);
% sippi_plot_prior(prior,m);
%
%
%
% % OPTIONAL OPTIONS
%
%   prior{1}.options.n_cond [int]: number of conditional points (def=5)
%   prior{1}.options.n_max_ite [int]: number of maximum iterations through the TI for matching patterns (def=200)
% 
%   prior{1}.options.plot    [int]: [0]:none, [1]:plot cond, [2]:storing movie (def=0)
%   prior{1}.options.verbose [int]: [0] no infor to screen, [1]:some info (def=1)
% 
% 
%
% TMH/2014
%
% See also: sippi_prior_init, sippi_prior
%
function [m_propose,prior]=sippi_prior_dsim(prior,m_current,ip);

if ~exist('mps_dsim.m','file')
    disp(sptinf('%s: mps_dsim.m is not in the path!!',mfilename));
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
if ~isfield(prior{ip}.options,'precalc_dist_full');
    prior{ip}.options.precalc_dist_full=0;
    try 
        if prod(prior{1}.dim)<=(80*80)
            prior{ip}.options.precalc_dist_full=1;
        end
    end
end

if nargin == 1
    if prior{ip}.ndim==1;
        SIM_data=NaN.*ones(1,length(prior{ip}.x));
    elseif prior{ip}.ndim==2;
        SIM_data=NaN.*ones(length(prior{ip}.y),length(prior{ip}.x));
    else
        SIM_data=NaN;
        disp(sprintf('%s: 3D not supported for DSIM yet',mfilename))
    end
else
    % PERTURB
    SIM_data=m_current{ip};
    N=prod(size(SIM_data));
    i_resample=randomsample(N,ceil(prior{1}.seq_gibbs.step*N));
    SIM_data(i_resample)=NaN;
end
[m_propose{ip},prior{ip}.options]=mps_dsim(prior{ip}.ti,SIM_data,prior{ip}.options);

%m_propose{ip}=normcdf(prior{ip}.randn,0,1)*(prior{ip}.max-prior{ip}.min)+prior{ip}.min;
    
